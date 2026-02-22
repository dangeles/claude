#!/usr/bin/env python3
"""
Claude Code Configuration Sync Tool

Bidirectionally syncs configuration between ~/.claude and this repository.

Usage:
    ./sync-config.py pull              # User-wide: ~/.claude -> repo
    ./sync-config.py push              # User-wide: repo -> ~/.claude
    ./sync-config.py pull-projects     # Project-specific configs -> repo
    ./sync-config.py push-projects     # Project-specific configs from repo
    ./sync-config.py status            # Show differences
    ./sync-config.py plan              # Create planning entry
    ./sync-config.py --dry-run pull    # Preview changes
"""

import argparse
import hashlib
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml


class Colors:
    """ANSI color codes for terminal output"""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class ConfigSync:
    def __init__(self, config_path: str = "sync.config.yaml", dry_run: bool = False,
                 verbose: bool = False, auto_yes: bool = False, delete_mode: bool = False):
        self.config_path = Path(config_path)
        self.dry_run = dry_run
        self.verbose = verbose
        self.auto_yes = auto_yes
        self.delete_mode = delete_mode
        self.repo_root = self.config_path.parent.resolve()

        # Load configuration
        with open(self.config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        # Expand paths
        self.source_dir = Path(self.config['source_dir']).expanduser()
        self.target_dir = (self.repo_root / self.config['target_dir']).resolve()

        # Lazy: resolve project config path but do NOT create directory yet
        project_cfg = self.config.get('project_configs', {})
        project_target = project_cfg.get('target_dir', '')
        if project_target:
            self.project_config_target = (self.repo_root / project_target).resolve()
        else:
            self.project_config_target = None

        # Create main target directory if needed (but NOT project_config_target)
        if not self.dry_run:
            self.target_dir.mkdir(parents=True, exist_ok=True)

    def _ensure_project_config_target(self):
        """Create project_config_target directory on first use (lazy initialization)."""
        if self.project_config_target and not self.dry_run:
            self.project_config_target.mkdir(parents=True, exist_ok=True)

    def print_info(self, message: str):
        """Print info message"""
        print(f"{Colors.OKBLUE}[INFO]{Colors.ENDC} {message}")

    def print_success(self, message: str):
        """Print success message"""
        print(f"{Colors.OKGREEN}[SUCCESS]{Colors.ENDC} {message}")

    def print_warning(self, message: str):
        """Print warning message"""
        print(f"{Colors.WARNING}[WARNING]{Colors.ENDC} {message}")

    def print_error(self, message: str):
        """Print error message"""
        print(f"{Colors.FAIL}[ERROR]{Colors.ENDC} {message}")

    def is_excluded(self, path: Path, relative_to: Path) -> bool:
        """Check if path matches any exclusion pattern.

        Exclusion patterns are matched against the relative path from relative_to.
        A pattern like 'debug/' or 'debug' excludes the directory and all contents.
        A pattern like 'settings.local.json' excludes that exact file.
        """
        try:
            rel_path = path.relative_to(relative_to)
        except ValueError:
            return False

        rel_str = str(rel_path)

        for exclusion in self.config['exclusions']:
            # Normalize: strip trailing slash for uniform matching
            exclusion_normalized = exclusion.rstrip('/')

            # Match: exact name, or path starts with exclusion as a directory prefix
            if rel_str == exclusion_normalized or rel_str.startswith(exclusion_normalized + '/'):
                return True

        return False

    def normalize_path_for_comparison(self, rel_path: Path) -> str:
        """Normalize a relative path for set comparison.

        On macOS (case-insensitive APFS), lowercases the path.
        On Linux (case-sensitive), returns as-is.
        """
        s = str(rel_path)
        if sys.platform == 'darwin':
            return s.lower()
        return s

    def walk_target_no_symlinks(self, target_dir: Path) -> List[Path]:
        """Walk target directory yielding regular files, never following symlinks.

        Skips symlinks entirely (both symlinked files and symlinked directories).
        Logs a warning for each symlink encountered.

        Returns: List of Path objects (absolute) that are regular files.
        """
        results = []
        if not target_dir.exists():
            return results

        for item in target_dir.rglob('*'):
            if item.is_symlink():
                self.print_warning(f"Skipping symlink: {item} -> {os.readlink(item)}")
                continue
            if item.is_file():
                results.append(item)
        return results

    def find_orphans(self, sync_rule_path: str) -> List[Tuple[Path, Path]]:
        """Find files in target (e.g. ~/.claude/<path>/) that have no counterpart in source (repo).

        Only operates on directory-type sync rules.

        Args:
            sync_rule_path: The path from the sync rule (e.g. "skills/")

        Returns:
            List of (absolute_target_path, relative_path) tuples for orphaned files.
            relative_path is in original case for display/deletion.
        """
        source_base = self.target_dir / sync_rule_path    # repo side: claude-config/<path>
        target_base = self.source_dir / sync_rule_path    # live side: ~/.claude/<path>

        if not source_base.is_dir() or not target_base.is_dir():
            return []

        # Build source file set (normalized for comparison)
        source_files_norm = set()
        for item in source_base.rglob('*'):
            if item.is_file() and not self.is_excluded(item, source_base):
                rel = item.relative_to(source_base)
                source_files_norm.add(self.normalize_path_for_comparison(rel))

        # Build target file set (walk without following symlinks)
        target_files = self.walk_target_no_symlinks(target_base)
        orphans = []
        for item in target_files:
            if self.is_excluded(item, target_base):
                continue
            rel = item.relative_to(target_base)
            norm = self.normalize_path_for_comparison(rel)
            if norm not in source_files_norm:
                orphans.append((item, rel))

        return orphans

    def delete_orphans(self, orphans: List[Tuple[Path, Path]], target_base: Path,
                       total_target_files: int) -> int:
        """Delete orphaned files with safety checks.

        Args:
            orphans: List of (absolute_path, relative_path) tuples from find_orphans
            target_base: The resolved base directory (e.g. ~/.claude/skills/)
            total_target_files: Total number of non-excluded files in target_base

        Returns:
            Number of files actually deleted (0 if dry_run)
        """
        if not orphans:
            return 0

        # Safety threshold (RISK-09): abort if deletions > 50% of target files
        if total_target_files > 0:
            ratio = len(orphans) / total_target_files
            if ratio > 0.5:
                self.print_error(
                    f"Safety threshold exceeded: {len(orphans)} of {total_target_files} files "
                    f"({ratio:.0%}) would be deleted. This usually indicates a misconfigured "
                    f"source path. Aborting. Use --dry-run to inspect."
                )
                sys.exit(1)

        # Print deletion manifest (always, even in --yes mode, for audit trail)
        print(f"\n{Colors.WARNING}Files to delete ({len(orphans)}):{Colors.ENDC}")
        for abs_path, rel_path in orphans:
            print(f"  {Colors.FAIL}-{Colors.ENDC} {rel_path}")

        # Dry run: stop here
        if self.dry_run:
            self.print_info(f"Would delete {len(orphans)} orphaned file(s)")
            return 0

        # Interactive confirmation (unless --yes)
        if not self.auto_yes:
            response = input(
                f"\n{Colors.BOLD}Delete {len(orphans)} orphaned file(s)? [y/N]:{Colors.ENDC} "
            ).strip().lower()
            if response != 'y':
                self.print_info("Deletion cancelled")
                return 0

        # Execute deletions
        deleted_count = 0
        resolved_base = target_base.resolve()

        for abs_path, rel_path in orphans:
            # Safety: verify path is still within target_base after resolution
            try:
                resolved = abs_path.resolve()
                if not resolved.is_relative_to(resolved_base):
                    self.print_warning(
                        f"Skipping {rel_path}: resolves outside target directory"
                    )
                    continue
            except OSError:
                self.print_warning(f"Skipping {rel_path}: cannot resolve path")
                continue

            # Safety: never delete symlinks (double-check; walk_target should have skipped them)
            if abs_path.is_symlink():
                self.print_warning(f"Skipping symlink: {rel_path}")
                continue

            try:
                abs_path.unlink()
                self.print_info(f"Deleted: {rel_path}")
                deleted_count += 1
            except OSError as e:
                self.print_error(f"Failed to delete {rel_path}: {e}")

        # Clean up empty directories (bottom-up)
        self._cleanup_empty_dirs(target_base)

        return deleted_count

    def _cleanup_empty_dirs(self, base_dir: Path):
        """Remove empty directories under base_dir, bottom-up. Never removes base_dir itself."""
        if not base_dir.is_dir():
            return
        # Walk bottom-up: sorted by depth (deepest first)
        dirs = sorted(
            [d for d in base_dir.rglob('*') if d.is_dir() and not d.is_symlink()],
            key=lambda p: len(p.parts),
            reverse=True
        )
        for d in dirs:
            if d == base_dir:
                continue
            try:
                d.rmdir()  # Only succeeds if empty
                if self.verbose:
                    self.print_info(f"Removed empty directory: {d.relative_to(base_dir)}")
            except OSError:
                pass  # Directory not empty, skip

    def validate_source_paths(self) -> bool:
        """Pre-flight check: verify all sync rule source directories exist.

        Returns True if all valid, False if any missing.
        When used with --delete, a missing source would cause all target files
        to be classified as orphans, so this is a critical safety gate.
        """
        all_valid = True
        for rule in self.config['sync_rules']['always']:
            path = rule['path']
            source = self.target_dir / path  # repo side
            if not source.exists():
                self.print_error(
                    f"Source path does not exist: {source} (sync rule: {path}). "
                    f"Cannot run --delete with missing sources."
                )
                all_valid = False
        return all_valid

    def compute_checksum(self, file_path: Path) -> str:
        """Compute SHA256 checksum of file"""
        sha256 = hashlib.sha256()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(8192), b''):
                sha256.update(chunk)
        return sha256.hexdigest()

    def files_differ(self, file1: Path, file2: Path) -> bool:
        """Check if two files have different content"""
        if not file1.exists() or not file2.exists():
            return True
        return self.compute_checksum(file1) != self.compute_checksum(file2)

    def show_diff(self, file1: Path, file2: Path, label1: str = "source", label2: str = "target"):
        """Show diff between two files using git diff or diff"""
        try:
            result = subprocess.run(
                ['git', 'diff', '--no-index', '--color=always', str(file2), str(file1)],
                capture_output=True,
                text=True
            )
            if result.stdout:
                print(result.stdout)
            elif result.stderr:
                # Git diff returns non-zero when files differ, check stderr
                if "diff --git" in result.stderr:
                    print(result.stderr)
        except FileNotFoundError:
            # Git not available, try basic diff
            try:
                result = subprocess.run(
                    ['diff', '-u', str(file2), str(file1)],
                    capture_output=True,
                    text=True
                )
                if result.stdout:
                    print(result.stdout)
            except FileNotFoundError:
                self.print_warning("diff tool not available, cannot show diff")

    def handle_conflict(self, source: Path, target: Path) -> str:
        """Handle file conflict with interactive prompt

        Returns: 'source', 'target', 'skip', or 'abort'
        """
        if self.auto_yes:
            self.print_info(f"[AUTO] Overwriting {target} with {source} (--yes mode)")
            return 'source'

        print(f"\n{Colors.WARNING}Conflict detected:{Colors.ENDC}")
        print(f"  Source: {source}")
        if source.exists():
            stat = source.stat()
            print(f"    Modified: {datetime.fromtimestamp(stat.st_mtime)}, Size: {stat.st_size} bytes")
        else:
            print(f"    {Colors.FAIL}Does not exist{Colors.ENDC}")

        print(f"  Target: {target}")
        if target.exists():
            stat = target.stat()
            print(f"    Modified: {datetime.fromtimestamp(stat.st_mtime)}, Size: {stat.st_size} bytes")
        else:
            print(f"    {Colors.FAIL}Does not exist{Colors.ENDC}")

        # Show diff if both files exist
        if source.exists() and target.exists() and self.config['conflict_resolution']['show_diff']:
            print(f"\n{Colors.BOLD}Diff preview:{Colors.ENDC}")
            self.show_diff(source, target)

        # Prompt for resolution
        print(f"\n{Colors.BOLD}Choose resolution:{Colors.ENDC}")
        print("[1] Use source (overwrite target)")
        print("[2] Use target (keep current)")
        print("[3] Show full diff")
        print("[4] Skip this file")
        if self.config['conflict_resolution']['allow_abort']:
            print("[5] Abort sync")

        while True:
            choice = input(f"\n{Colors.BOLD}Choice:{Colors.ENDC} ").strip()

            if choice == '1':
                return 'source'
            elif choice == '2':
                return 'target'
            elif choice == '3':
                if source.exists() and target.exists():
                    self.show_diff(source, target)
                else:
                    self.print_warning("Cannot show diff - one or both files missing")
            elif choice == '4':
                return 'skip'
            elif choice == '5' and self.config['conflict_resolution']['allow_abort']:
                return 'abort'
            else:
                print(f"{Colors.FAIL}Invalid choice. Please try again.{Colors.ENDC}")

    def create_backup(self, file_path: Path) -> Optional[Path]:
        """Create backup of file before overwriting"""
        if not file_path.exists():
            return None

        backup_dir = self.repo_root / self.config['backup']['location']
        backup_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_name = f"{file_path.name}.{timestamp}.bak"
        backup_path = backup_dir / backup_name

        shutil.copy2(file_path, backup_path)
        self.print_info(f"Created backup: {backup_path}")

        # Cleanup old backups
        self.cleanup_old_backups(backup_dir, file_path.name)

        return backup_path

    def cleanup_old_backups(self, backup_dir: Path, file_name: str):
        """Keep only max_backups most recent backups"""
        max_backups = self.config['backup']['max_backups']

        # Find all backups for this file
        backups = sorted(
            [f for f in backup_dir.glob(f"{file_name}.*.bak")],
            key=lambda p: p.stat().st_mtime,
            reverse=True
        )

        # Remove old backups
        for backup in backups[max_backups:]:
            backup.unlink()
            if self.verbose:
                self.print_info(f"Removed old backup: {backup}")

    def copy_file(self, source: Path, target: Path, check_conflict: bool = True) -> bool:
        """Copy file from source to target with conflict handling

        Returns: True if file was copied, False if skipped
        """
        # Check if target exists and differs
        if target.exists() and check_conflict:
            if not self.files_differ(source, target):
                if self.verbose:
                    self.print_info(f"Skipping {target.name} (identical)")
                return False

            # Handle conflict
            resolution = self.handle_conflict(source, target)

            if resolution == 'abort':
                self.print_error("Sync aborted by user")
                sys.exit(1)
            elif resolution == 'skip':
                self.print_info(f"Skipped: {target.name}")
                return False
            elif resolution == 'target':
                self.print_info(f"Kept target: {target.name}")
                return False
            # resolution == 'source' -> continue to copy

        # Create backup if enabled
        if self.config['backup']['enabled'] and target.exists():
            self.create_backup(target)

        # Copy file
        if self.dry_run:
            self.print_info(f"Would copy: {source} -> {target}")
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            self.print_success(f"Copied: {source.name}")

        return True

    def copy_tree(self, source: Path, target: Path, check_conflict: bool = True):
        """Recursively copy directory tree with exclusions"""
        if not source.exists():
            self.print_warning(f"Source does not exist: {source}")
            return

        for item in source.rglob('*'):
            if item.is_file():
                # Check exclusions
                if self.is_excluded(item, source):
                    if self.verbose:
                        self.print_info(f"Excluded: {item.relative_to(source)}")
                    continue

                # Compute target path
                rel_path = item.relative_to(source)
                target_file = target / rel_path

                # Copy file
                self.copy_file(item, target_file, check_conflict)

    def pull_config(self):
        """Pull configuration from ~/.claude to repo"""
        self.print_info(f"Pulling configuration from {self.source_dir} to {self.target_dir}")

        if not self.source_dir.exists():
            self.print_error(f"Source directory does not exist: {self.source_dir}")
            sys.exit(1)

        # Sync each configured item
        for rule in self.config['sync_rules']['always']:
            path = rule['path']
            source = self.source_dir / path
            target = self.target_dir / path

            if not source.exists():
                self.print_warning(f"Source not found: {source}")
                continue

            self.print_info(f"Syncing: {path}")

            if source.is_dir():
                self.copy_tree(source, target, check_conflict=True)
            else:
                self.copy_file(source, target, check_conflict=True)

        self.print_success("Pull complete!")

    def push_config(self):
        """Push configuration from repo to ~/.claude"""
        self.print_warning(f"About to modify system Claude Code configuration at: {self.source_dir}")

        if self.delete_mode:
            self.print_warning("--delete is enabled: orphaned files will be removed from target")

        if not self.dry_run and not self.auto_yes:
            response = input(f"{Colors.BOLD}Continue? [y/N]:{Colors.ENDC} ").strip().lower()
            if response != 'y':
                self.print_info("Push cancelled")
                return

        # Pre-flight validation for --delete
        if self.delete_mode:
            if not self.validate_source_paths():
                self.print_error("Aborting: fix source paths before using --delete")
                sys.exit(1)

        self.print_info(f"Pushing configuration from {self.target_dir} to {self.source_dir}")

        # --- Copy phase (existing logic) ---
        for rule in self.config['sync_rules']['always']:
            path = rule['path']
            source = self.target_dir / path
            target = self.source_dir / path

            if not source.exists():
                self.print_warning(f"Source not found: {source}")
                continue

            self.print_info(f"Syncing: {path}")

            if source.is_dir():
                self.copy_tree(source, target, check_conflict=True)
            else:
                self.copy_file(source, target, check_conflict=True)

        # --- Deletion phase (only if --delete) ---
        if self.delete_mode:
            total_deleted = 0
            for rule in self.config['sync_rules']['always']:
                path = rule['path']
                source_base = self.target_dir / path
                target_base = self.source_dir / path

                # Only process directory-type rules for orphan detection
                if not source_base.is_dir():
                    continue

                self.print_info(f"Scanning for orphans: {path}")
                orphans = self.find_orphans(path)

                if orphans:
                    # Count total non-excluded files in target for threshold calculation
                    target_files = self.walk_target_no_symlinks(target_base)
                    non_excluded = [f for f in target_files if not self.is_excluded(f, target_base)]
                    total_deleted += self.delete_orphans(orphans, target_base, len(non_excluded))
                else:
                    if self.verbose:
                        self.print_info(f"No orphans found in {path}")

            if total_deleted > 0:
                self.print_success(f"Deleted {total_deleted} orphaned file(s)")

        self.print_success("Push complete!")

    def discover_projects(self) -> List[Tuple[str, Path]]:
        """Discover projects with Claude Code configuration

        Returns: List of (project_name, project_path) tuples
        """
        projects = []

        for search_path in self.config['project_configs']['search_paths']:
            search_path = Path(search_path).expanduser()

            # Handle wildcards
            if '*' in str(search_path):
                parent = Path(str(search_path).split('*')[0])
                pattern = str(search_path).split('*')[1] if '*' in str(search_path).split('*')[1] else ''

                if not parent.exists():
                    continue

                for project_dir in parent.iterdir():
                    if not project_dir.is_dir():
                        continue

                    # Check if project has Claude Code config
                    has_config = False
                    for config_file in self.config['project_configs']['config_files']:
                        if (project_dir / config_file).exists():
                            has_config = True
                            break

                    if has_config:
                        projects.append((project_dir.name, project_dir))
            else:
                # Exact path
                if search_path.exists():
                    has_config = False
                    for config_file in self.config['project_configs']['config_files']:
                        if (search_path / config_file).exists():
                            has_config = True
                            break

                    if has_config:
                        projects.append((search_path.name, search_path))

        return projects

    def pull_project_configs(self, project_name: Optional[str] = None):
        """Pull project-specific configurations to repo"""
        self._ensure_project_config_target()

        if self.project_config_target is None:
            self.print_warning("project_configs.target_dir is not configured. Skipping.")
            return

        self.print_info("Discovering projects with Claude Code configuration...")

        projects = self.discover_projects()

        if not projects:
            self.print_warning("No projects found with Claude Code configuration")
            return

        self.print_info(f"Found {len(projects)} project(s):")
        for name, path in projects:
            print(f"  - {name} ({path})")

        # Filter by project name if specified
        if project_name:
            projects = [(n, p) for n, p in projects if n == project_name]
            if not projects:
                self.print_error(f"Project not found: {project_name}")
                return

        # Sync each project
        for name, path in projects:
            self.print_info(f"\nSyncing project: {name}")

            project_target = self.project_config_target / name
            project_target.mkdir(parents=True, exist_ok=True)

            for config_file in self.config['project_configs']['config_files']:
                source = path / config_file

                if not source.exists():
                    continue

                if source.is_dir():
                    target = project_target / config_file
                    self.copy_tree(source, target, check_conflict=True)
                else:
                    target = project_target / config_file
                    self.copy_file(source, target, check_conflict=True)

        self.print_success("Project config pull complete!")

    def push_project_configs(self, project_name: Optional[str] = None):
        """Push project-specific configurations from repo"""
        self._ensure_project_config_target()

        if self.project_config_target is None:
            self.print_warning("project_configs.target_dir is not configured. Skipping.")
            return

        self.print_info("Finding project configurations in repository...")

        if not self.project_config_target.exists():
            self.print_warning("No project configurations found in repository")
            return

        # Find all project directories
        project_dirs = [d for d in self.project_config_target.iterdir() if d.is_dir()]

        if not project_dirs:
            self.print_warning("No project configurations found in repository")
            return

        # Filter by project name if specified
        if project_name:
            project_dirs = [d for d in project_dirs if d.name == project_name]
            if not project_dirs:
                self.print_error(f"Project configuration not found: {project_name}")
                return

        self.print_info(f"Found {len(project_dirs)} project configuration(s)")

        # Sync each project
        for project_dir in project_dirs:
            name = project_dir.name
            self.print_info(f"\nPushing configuration for project: {name}")

            # Try to find the project in search paths
            project_path = None
            for search_path in self.config['project_configs']['search_paths']:
                search_path = Path(search_path).expanduser()

                if '*' in str(search_path):
                    parent = Path(str(search_path).split('*')[0])
                    potential_path = parent / name
                    if potential_path.exists():
                        project_path = potential_path
                        break

            if not project_path:
                self.print_warning(f"Project directory not found: {name}")
                self.print_info("  Skipping (project may not exist on this machine)")
                continue

            self.print_info(f"  Target: {project_path}")

            # Copy config files
            for config_file in self.config['project_configs']['config_files']:
                source = project_dir / config_file

                if not source.exists():
                    continue

                target = project_path / config_file

                if source.is_dir():
                    self.copy_tree(source, target, check_conflict=True)
                else:
                    self.copy_file(source, target, check_conflict=True)

        # --- Deletion phase for project configs (only if --delete) ---
        if self.delete_mode:
            total_deleted = 0
            for project_dir in project_dirs:
                name = project_dir.name

                # Find the project's live path
                project_path = None
                for search_path in self.config['project_configs']['search_paths']:
                    sp = Path(search_path).expanduser()
                    if '*' in str(sp):
                        parent = Path(str(sp).split('*')[0])
                        potential_path = parent / name
                        if potential_path.exists():
                            project_path = potential_path
                            break

                if not project_path:
                    continue

                # Check each config directory for orphans
                for config_file in self.config['project_configs']['config_files']:
                    source_dir = project_dir / config_file
                    target_dir = project_path / config_file

                    if not source_dir.is_dir() or not target_dir.is_dir():
                        continue

                    # Build source file set
                    source_files_norm = set()
                    for item in source_dir.rglob('*'):
                        if item.is_file() and not self.is_excluded(item, source_dir):
                            rel = item.relative_to(source_dir)
                            source_files_norm.add(self.normalize_path_for_comparison(rel))

                    # Build target file set
                    target_files = self.walk_target_no_symlinks(target_dir)
                    orphans = []
                    for item in target_files:
                        if self.is_excluded(item, target_dir):
                            continue
                        rel = item.relative_to(target_dir)
                        norm = self.normalize_path_for_comparison(rel)
                        if norm not in source_files_norm:
                            orphans.append((item, rel))

                    if orphans:
                        non_excluded = [f for f in target_files if not self.is_excluded(f, target_dir)]
                        self.print_info(f"Scanning for orphans in project {name}/{config_file}")
                        total_deleted += self.delete_orphans(orphans, target_dir, len(non_excluded))

            if total_deleted > 0:
                self.print_success(f"Deleted {total_deleted} orphaned project config file(s)")

        self.print_success("Project config push complete!")

    def show_status(self):
        """Show differences between source and target"""
        self.print_info("Comparing configuration...")

        print(f"\n{Colors.BOLD}User-wide configuration:{Colors.ENDC}")
        print(f"  Source: {self.source_dir}")
        print(f"  Target: {self.target_dir}")

        changes = []

        for rule in self.config['sync_rules']['always']:
            path = rule['path']
            source = self.source_dir / path
            target = self.target_dir / path

            if not source.exists() and not target.exists():
                continue

            if not source.exists():
                changes.append((path, "deleted in source"))
            elif not target.exists():
                changes.append((path, "new in source"))
            elif source.is_file() and target.is_file():
                if self.files_differ(source, target):
                    changes.append((path, "modified"))
            elif source.is_dir() and target.is_dir():
                # Check for differences in directory
                source_files = {f.relative_to(source) for f in source.rglob('*') if f.is_file() and not self.is_excluded(f, source)}
                target_files = {f.relative_to(target) for f in target.rglob('*') if f.is_file()}

                new_files = source_files - target_files
                deleted_files = target_files - source_files
                common_files = source_files & target_files

                modified_files = []
                for rel_path in common_files:
                    if self.files_differ(source / rel_path, target / rel_path):
                        modified_files.append(rel_path)

                if new_files or deleted_files or modified_files:
                    detail = []
                    if new_files:
                        detail.append(f"{len(new_files)} new")
                    if deleted_files:
                        detail.append(f"{len(deleted_files)} deleted")
                    if modified_files:
                        detail.append(f"{len(modified_files)} modified")
                    changes.append((path, ", ".join(detail)))

        if changes:
            print(f"\n{Colors.WARNING}Changes detected:{Colors.ENDC}")
            for path, status in changes:
                print(f"  {Colors.WARNING}*{Colors.ENDC} {path}: {status}")
        else:
            print(f"\n{Colors.OKGREEN}No changes detected{Colors.ENDC}")

        # --- Orphan detection (files in ~/.claude/ but not in repo) ---
        all_orphans = []
        for rule in self.config['sync_rules']['always']:
            path = rule['path']
            source_base = self.target_dir / path
            target_base = self.source_dir / path

            if not source_base.is_dir() or not target_base.is_dir():
                continue

            orphans = self.find_orphans(path)
            for abs_path, rel_path in orphans:
                all_orphans.append((path, rel_path))

        if all_orphans:
            print(f"\n{Colors.WARNING}Orphaned files ({len(all_orphans)}):{Colors.ENDC}")
            print(f"  (present in {self.source_dir} but not in {self.target_dir})")
            for rule_path, rel_path in all_orphans:
                print(f"  {Colors.FAIL}-{Colors.ENDC} {rule_path}{rel_path}")
            print(f"\n  Use {Colors.BOLD}push --delete{Colors.ENDC} to remove orphaned files")
        else:
            if self.verbose:
                print(f"\n{Colors.OKGREEN}No orphaned files detected{Colors.ENDC}")

        # Project configurations
        print(f"\n{Colors.BOLD}Project configurations:{Colors.ENDC}")
        projects = self.discover_projects()
        if projects:
            print(f"  Found {len(projects)} project(s) with configuration:")
            for name, path in projects:
                print(f"    {Colors.OKCYAN}*{Colors.ENDC} {name} ({path})")
        else:
            print("  No projects found")

    def detect_hostname(self) -> str:
        """Detect machine hostname"""
        return platform.node().split('.')[0]  # Remove domain if present

    def create_plan_entry(self, title: str, status: str = "Planned", open_editor: bool = False):
        """Create new planning journal entry"""
        hostname = self.detect_hostname()
        date = datetime.now().strftime("%Y-%m-%d")

        # Create machine directory
        machine_dir = self.repo_root / "planning" / hostname
        machine_dir.mkdir(parents=True, exist_ok=True)

        # Sanitize title for filename
        filename_title = title.lower().replace(' ', '-').replace('/', '-')
        filename_title = ''.join(c for c in filename_title if c.isalnum() or c == '-')

        filename = f"{date}-{filename_title}.md"
        filepath = machine_dir / filename

        # Load template
        template_path = self.repo_root / "planning" / ".template.md"
        with open(template_path, 'r') as f:
            template = f.read()

        # Fill in template
        content = template.replace("YYYY-MM-DD", date)
        content = content.replace("[auto-detected hostname]", hostname)
        content = content.replace("[Planned / In Progress / Success / Partial / Failed]", status)
        content = content.replace("[Brief Description of Config Change]", title)

        # Write file
        if self.dry_run:
            self.print_info(f"Would create: {filepath}")
        else:
            with open(filepath, 'w') as f:
                f.write(content)
            self.print_success(f"Created planning entry: {filepath}")

            # Open in editor only if --edit was passed
            if open_editor:
                editor = os.environ.get('EDITOR', 'nano')
                try:
                    subprocess.run([editor, str(filepath)])
                except FileNotFoundError:
                    self.print_info(f"Edit file at: {filepath}")
            else:
                self.print_info(f"Created: {filepath}. Use --edit to open in editor.")

        return filepath

    def list_plans(self, machine: Optional[str] = None):
        """List planning journal entries"""
        planning_dir = self.repo_root / "planning"

        if not planning_dir.exists():
            self.print_warning("No planning entries found")
            return

        # Collect all entries
        entries = []

        for machine_dir in planning_dir.iterdir():
            if not machine_dir.is_dir() or machine_dir.name.startswith('.'):
                continue

            if machine and machine_dir.name != machine:
                continue

            for entry in machine_dir.glob("*.md"):
                entries.append((machine_dir.name, entry))

        if not entries:
            self.print_warning("No planning entries found")
            return

        # Sort by date (newest first)
        entries.sort(key=lambda e: e[1].name, reverse=True)

        print(f"\n{Colors.BOLD}Planning Journal Entries:{Colors.ENDC}")

        current_machine = None
        for machine, entry in entries:
            if machine != current_machine:
                print(f"\n{Colors.OKCYAN}{machine}:{Colors.ENDC}")
                current_machine = machine

            # Extract title from filename (date-title.md)
            parts = entry.stem.split('-', 3)
            if len(parts) >= 4:
                date = '-'.join(parts[:3])
                title = parts[3].replace('-', ' ').title()
            else:
                date = entry.stem
                title = "Untitled"

            print(f"  {date}: {title}")
            print(f"    {Colors.OKBLUE}{entry}{Colors.ENDC}")


def main():
    parser = argparse.ArgumentParser(
        description="Claude Code Configuration Sync Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ./sync-config.py pull                    # Pull user-wide config
  ./sync-config.py push                    # Push user-wide config
  ./sync-config.py push --yes              # Push without prompts
  ./sync-config.py push --yes --delete     # Push and remove orphaned files
  ./sync-config.py pull-projects           # Pull all project configs
  ./sync-config.py push-projects bioreactor  # Push specific project config
  ./sync-config.py status                  # Show differences
  ./sync-config.py plan --title "Enable new plugin"
  ./sync-config.py plan --title "Enable new plugin" --edit
  ./sync-config.py plan --list             # List all plans
        """
    )

    parser.add_argument('command',
                        choices=['pull', 'push', 'pull-projects', 'push-projects', 'status', 'plan'],
                        help="Command to execute")
    parser.add_argument('project', nargs='?', help="Project name (for push-projects/pull-projects)")
    parser.add_argument('--dry-run', action='store_true', help="Preview changes without executing")
    parser.add_argument('--verbose', '-v', action='store_true', help="Verbose output")
    parser.add_argument('--config', default='sync.config.yaml', help="Config file path")
    parser.add_argument('--no-backup', action='store_true', help="Disable backups")
    parser.add_argument('--yes', '-y', action='store_true',
                        help="Skip all interactive prompts; source always wins conflicts")
    parser.add_argument('--delete', action='store_true',
                        help="Remove orphaned files in target that have no source counterpart "
                             "(only within sync_rules paths)")

    # Plan command options
    parser.add_argument('--title', help="Title for new planning entry")
    parser.add_argument('--status', help="Status for planning entry")
    parser.add_argument('--list', dest='list_plans', action='store_true', help="List planning entries")
    parser.add_argument('--machine', help="Filter by machine hostname")
    parser.add_argument('--edit', action='store_true',
                        help="Open editor after creating plan entry (default: just create file)")

    args = parser.parse_args()

    # Create sync instance
    sync = ConfigSync(args.config, dry_run=args.dry_run, verbose=args.verbose,
                      auto_yes=args.yes, delete_mode=args.delete)

    # Disable backups if requested
    if args.no_backup:
        sync.config['backup']['enabled'] = False

    # Execute command
    try:
        if args.command == 'pull':
            sync.pull_config()
        elif args.command == 'push':
            sync.push_config()
        elif args.command == 'pull-projects':
            sync.pull_project_configs(args.project)
        elif args.command == 'push-projects':
            sync.push_project_configs(args.project)
        elif args.command == 'status':
            sync.show_status()
        elif args.command == 'plan':
            if args.list_plans:
                sync.list_plans(args.machine)
            else:
                if not args.title:
                    parser.error("--title required for creating planning entry")
                status = args.status or "Planned"
                sync.create_plan_entry(args.title, status, open_editor=args.edit)

    except KeyboardInterrupt:
        print(f"\n{Colors.WARNING}Interrupted by user{Colors.ENDC}")
        sys.exit(1)
    except Exception as e:
        print(f"\n{Colors.FAIL}Error: {e}{Colors.ENDC}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
