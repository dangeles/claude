#!/usr/bin/env python3
"""Handoff validation script for programming-pm workflow.

Validates handoffs against handoff-schema.md specification (version 1.1).
Supports 8 handoff types across all phase boundaries.

Usage:
    python3 validate-handoff.py <handoff_path> <handoff_type>

Example:
    python3 validate-handoff.py /tmp/session/handoffs/phase1-requirements-handoff.yaml requirements_handoff

Exit codes:
    0: Validation successful
    1: Validation failed (errors found)
    2: Invalid usage
"""

import sys
import yaml
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional


class HandoffValidator:
    """Base validator with common validation logic."""

    def __init__(self, handoff_data: Dict[str, Any], handoff_path: Path):
        self.data = handoff_data
        self.path = handoff_path
        self.errors: List[str] = []

    def validate_base_fields(self) -> None:
        """Validate base handoff fields present in all handoffs."""
        # Check handoff wrapper
        if 'handoff' not in self.data:
            self.errors.append("Missing root 'handoff' field")
            return

        handoff = self.data['handoff']

        # Required metadata
        self._check_field(handoff, 'version', str)
        if 'version' in handoff and handoff['version'] != '1.1':
            self.errors.append(f"Invalid version: {handoff['version']} (expected 1.1)")

        self._check_field(handoff, 'from_phase', int)
        self._check_field(handoff, 'to_phase', int)
        self._check_field(handoff, 'producer', str)
        self._check_field(handoff, 'consumer', str)
        self._check_field(handoff, 'timestamp', str)

        # Validate ISO8601 timestamp
        if 'timestamp' in handoff:
            try:
                datetime.fromisoformat(handoff['timestamp'].replace('Z', '+00:00'))
            except (ValueError, AttributeError):
                self.errors.append(f"Invalid ISO8601 timestamp: {handoff.get('timestamp')}")

        # Session context
        if 'session' in handoff:
            session = handoff['session']
            self._check_field(session, 'session_dir', str)
            if 'session_dir' in session and not session['session_dir'].startswith('/'):
                self.errors.append(f"session_dir must be absolute path: {session['session_dir']}")

            self._check_field(session, 'archival_guidelines_path', str)
            if 'archival_guidelines_path' in session and not session['archival_guidelines_path'].startswith('/'):
                self.errors.append(f"archival_guidelines_path must be absolute path")
        else:
            self.errors.append("Missing 'session' field")

        # Deliverable reference
        if 'deliverable' in handoff:
            deliverable = handoff['deliverable']
            self._check_field(deliverable, 'location', str)
            if 'location' in deliverable and not deliverable['location'].startswith('/'):
                self.errors.append(f"deliverable.location must be absolute path")

            self._check_field(deliverable, 'type', str)
            valid_types = ['specification', 'code', 'tests', 'documentation']
            if 'type' in deliverable and deliverable['type'] not in valid_types:
                self.errors.append(f"Invalid deliverable type: {deliverable['type']}")
        else:
            self.errors.append("Missing 'deliverable' field")

        # Context
        if 'context' in handoff:
            context = handoff['context']
            self._check_field(context, 'task_id', str)
            self._check_field(context, 'description', str)
            self._check_field(context, 'focus_areas', list)
            self._check_field(context, 'known_gaps', list)
        else:
            self.errors.append("Missing 'context' field")

        # Quality assessment
        if 'quality' in handoff:
            quality = handoff['quality']
            self._check_field(quality, 'status', str)
            if 'status' in quality and quality['status'] not in ['complete', 'partial']:
                self.errors.append(f"Invalid quality status: {quality['status']}")

            self._check_field(quality, 'confidence', str)
            if 'confidence' in quality and quality['confidence'] not in ['high', 'medium', 'low']:
                self.errors.append(f"Invalid confidence level: {quality['confidence']}")

            self._check_field(quality, 'notes', str)
        else:
            self.errors.append("Missing 'quality' field")

    def _check_field(self, obj: Dict[str, Any], field: str, expected_type: type) -> None:
        """Check if field exists and has correct type."""
        if field not in obj:
            self.errors.append(f"Missing required field: {field}")
        elif not isinstance(obj[field], expected_type):
            self.errors.append(f"Invalid type for {field}: expected {expected_type.__name__}, got {type(obj[field]).__name__}")

    def _check_optional_field(self, obj: Dict[str, Any], field: str, expected_type: type) -> None:
        """Check if optional field has correct type when present."""
        if field in obj and not isinstance(obj[field], expected_type):
            self.errors.append(f"Invalid type for {field}: expected {expected_type.__name__}, got {type(obj[field]).__name__}")

    def validate(self) -> List[str]:
        """Override in subclasses to add specific validation."""
        self.validate_base_fields()
        return self.errors


class SessionHandoffValidator(HandoffValidator):
    """Validator for Phase 0 -> Phase 1 (session_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Session-specific fields
        if 'session' in handoff:
            session = handoff['session']
            self._check_field(session, 'guidelines_found', bool)
            self._check_field(session, 'guidelines_source', str)

        # Archival guidelines
        if 'archival_guidelines' in handoff:
            guidelines = handoff['archival_guidelines']

            self._check_field(guidelines, 'code_directories', list)
            if 'code_directories' in guidelines:
                for i, dir_entry in enumerate(guidelines['code_directories']):
                    self._check_field(dir_entry, 'name', str)
                    self._check_field(dir_entry, 'purpose', str)

            if 'git_workflow' in guidelines:
                git = guidelines['git_workflow']
                self._check_field(git, 'commit_after_edit', bool)
                self._check_field(git, 'stage_specific_files', bool)
                self._check_field(git, 'no_destructive_ops', bool)
                self._check_field(git, 'conventional_commits', bool)

            if 'testing_conventions' in guidelines:
                testing = guidelines['testing_conventions']
                self._check_optional_field(testing, 'coverage_target', (int, float))
                self._check_optional_field(testing, 'test_directory', str)

            if 'documentation_conventions' in guidelines:
                docs = guidelines['documentation_conventions']
                self._check_field(docs, 'docstrings_required', bool)
                self._check_field(docs, 'type_hints_required', bool)
                self._check_field(docs, 'readme_updates', bool)

            if 'code_style' in guidelines:
                style = guidelines['code_style']
                self._check_field(style, 'linter', str)
                self._check_optional_field(style, 'formatter', str)
        else:
            self.errors.append("Missing 'archival_guidelines' field")

        return self.errors


class RequirementsHandoffValidator(HandoffValidator):
    """Validator for Phase 1 -> Phase 2 (requirements_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Requirements
        if 'requirements' in handoff:
            req = handoff['requirements']
            self._check_field(req, 'problem_statement', str)
            self._check_field(req, 'success_criteria', list)

            if 'scope' in req:
                scope = req['scope']
                self._check_field(scope, 'in_scope', list)
                self._check_field(scope, 'out_of_scope', list)
            else:
                self.errors.append("Missing 'requirements.scope' field")

            self._check_field(req, 'constraints', list)
            self._check_field(req, 'dependencies', list)
        else:
            self.errors.append("Missing 'requirements' field")

        # Stakeholders
        if 'stakeholders' in handoff:
            stakeholders = handoff['stakeholders']
            self._check_field(stakeholders, 'primary', str)
            self._check_field(stakeholders, 'consulted', list)
        else:
            self.errors.append("Missing 'stakeholders' field")

        # Approval
        if 'approval' in handoff:
            approval = handoff['approval']
            self._check_field(approval, 'approved_by', str)
            self._check_field(approval, 'approved_date', str)
            if 'approved_date' in approval:
                try:
                    datetime.fromisoformat(approval['approved_date'].replace('Z', '+00:00'))
                except (ValueError, AttributeError):
                    self.errors.append(f"Invalid ISO8601 approved_date")
            self._check_field(approval, 'conditions', list)
        else:
            self.errors.append("Missing 'approval' field")

        return self.errors


class PremortемHandoffValidator(HandoffValidator):
    """Validator for Phase 2 -> Phase 3 (premortem_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Risks
        if 'risks' in handoff:
            risks = handoff['risks']
            self._check_field(risks, 'identified', int)
            self._check_field(risks, 'critical', list)
            self._check_field(risks, 'high', list)
            self._check_field(risks, 'mitigated', list)
            self._check_field(risks, 'accepted', list)
        else:
            self.errors.append("Missing 'risks' field")

        # Risk summary
        if 'risk_summary' in handoff:
            for i, risk in enumerate(handoff['risk_summary']):
                self._check_field(risk, 'id', str)
                self._check_field(risk, 'description', str)
                self._check_field(risk, 'score', int)
                self._check_field(risk, 'disposition', str)
                if 'disposition' in risk and risk['disposition'] not in ['mitigate', 'accept', 'transfer', 'avoid']:
                    self.errors.append(f"Invalid disposition in risk {i}: {risk['disposition']}")
                if 'disposition' in risk and risk['disposition'] == 'mitigate':
                    self._check_field(risk, 'mitigation', str)
        else:
            self.errors.append("Missing 'risk_summary' field")

        # Architecture implications
        if 'architecture_implications' in handoff:
            for i, impl in enumerate(handoff['architecture_implications']):
                self._check_field(impl, 'risk_id', str)
                self._check_field(impl, 'implication', str)
        else:
            self.errors.append("Missing 'architecture_implications' field")

        return self.errors


class ArchitectureHandoffValidator(HandoffValidator):
    """Validator for Phase 3 -> Phase 4 (architecture_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Components
        if 'components' in handoff:
            for i, comp in enumerate(handoff['components']):
                self._check_field(comp, 'name', str)
                self._check_field(comp, 'responsibility', str)

                if 'interfaces' in comp:
                    interfaces = comp['interfaces']
                    self._check_field(interfaces, 'inputs', list)
                    self._check_field(interfaces, 'outputs', list)
                else:
                    self.errors.append(f"Missing 'interfaces' in component {i}")

                self._check_field(comp, 'dependencies', list)
                self._check_field(comp, 'estimated_effort', str)
        else:
            self.errors.append("Missing 'components' field")

        # Data flow
        if 'data_flow' in handoff:
            data_flow = handoff['data_flow']
            self._check_field(data_flow, 'description', str)
            self._check_optional_field(data_flow, 'diagram_location', str)
        else:
            self.errors.append("Missing 'data_flow' field")

        # Technology choices
        if 'technology_choices' in handoff:
            for i, tech in enumerate(handoff['technology_choices']):
                self._check_field(tech, 'category', str)
                self._check_field(tech, 'choice', str)
                self._check_field(tech, 'rationale', str)
        else:
            self.errors.append("Missing 'technology_choices' field")

        # Testing strategy
        if 'testing_strategy' in handoff:
            testing = handoff['testing_strategy']
            self._check_field(testing, 'unit', str)
            self._check_field(testing, 'integration', str)
            self._check_field(testing, 'coverage_target', (int, float))
        else:
            self.errors.append("Missing 'testing_strategy' field")

        # Implementation order
        if 'implementation_order' in handoff:
            for i, item in enumerate(handoff['implementation_order']):
                self._check_field(item, 'component', str)
                self._check_field(item, 'priority', int)
                self._check_field(item, 'dependencies', list)
        else:
            self.errors.append("Missing 'implementation_order' field")

        return self.errors


class MathHandoffValidator(HandoffValidator):
    """Validator for mathematician -> developer (math_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Algorithm
        if 'algorithm' in handoff:
            algo = handoff['algorithm']
            self._check_field(algo, 'name', str)
            self._check_field(algo, 'description', str)
            self._check_field(algo, 'pseudocode', str)
        else:
            self.errors.append("Missing 'algorithm' field")

        # Complexity analysis
        if 'complexity_analysis' in handoff:
            complexity = handoff['complexity_analysis']

            if 'time' in complexity:
                time = complexity['time']
                self._check_field(time, 'best_case', str)
                self._check_field(time, 'average_case', str)
                self._check_field(time, 'worst_case', str)
            else:
                self.errors.append("Missing 'complexity_analysis.time' field")

            if 'space' in complexity:
                space = complexity['space']
                self._check_field(space, 'auxiliary', str)
                self._check_field(space, 'total', str)
            else:
                self.errors.append("Missing 'complexity_analysis.space' field")

            self._check_field(complexity, 'analysis_notes', str)
        else:
            self.errors.append("Missing 'complexity_analysis' field")

        # Numerical stability
        if 'numerical_stability' in handoff:
            stability = handoff['numerical_stability']
            self._check_field(stability, 'stable', bool)
            self._check_field(stability, 'conditions', str)
            self._check_field(stability, 'precision_requirements', str)

            if 'failure_modes' in stability:
                for i, mode in enumerate(stability['failure_modes']):
                    self._check_field(mode, 'condition', str)
                    self._check_field(mode, 'symptom', str)
                    self._check_field(mode, 'mitigation', str)
            else:
                self.errors.append("Missing 'numerical_stability.failure_modes' field")
        else:
            self.errors.append("Missing 'numerical_stability' field")

        # Implementation guidance
        if 'implementation_guidance' in handoff:
            guidance = handoff['implementation_guidance']
            self._check_field(guidance, 'recommended_approach', str)

            if 'libraries' in guidance:
                for i, lib in enumerate(guidance['libraries']):
                    self._check_field(lib, 'name', str)
                    self._check_field(lib, 'usage', str)
            else:
                self.errors.append("Missing 'implementation_guidance.libraries' field")

            self._check_field(guidance, 'pitfalls', list)
        else:
            self.errors.append("Missing 'implementation_guidance' field")

        # Verification criteria
        if 'verification_criteria' in handoff:
            verification = handoff['verification_criteria']
            self._check_field(verification, 'invariants', list)

            if 'test_cases' in verification:
                for i, test in enumerate(verification['test_cases']):
                    self._check_field(test, 'name', str)
                    self._check_field(test, 'input', str)
                    self._check_field(test, 'expected', str)
            else:
                self.errors.append("Missing 'verification_criteria.test_cases' field")

            if 'edge_cases' in verification:
                for i, edge in enumerate(verification['edge_cases']):
                    self._check_field(edge, 'name', str)
                    self._check_field(edge, 'input', str)
                    self._check_field(edge, 'expected', str)
                    self._check_field(edge, 'note', str)
            else:
                self.errors.append("Missing 'verification_criteria.edge_cases' field")
        else:
            self.errors.append("Missing 'verification_criteria' field")

        return self.errors


class StatsHandoffValidator(HandoffValidator):
    """Validator for statistician -> developer (stats_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Method
        if 'method' in handoff:
            method = handoff['method']
            self._check_field(method, 'name', str)
            self._check_field(method, 'description', str)
            self._check_field(method, 'rationale', str)
        else:
            self.errors.append("Missing 'method' field")

        # Assumptions
        if 'assumptions' in handoff:
            assumptions = handoff['assumptions']
            self._check_field(assumptions, 'data_requirements', list)
            self._check_field(assumptions, 'distributional', list)

            if 'violations_impact' in assumptions:
                for i, impact in enumerate(assumptions['violations_impact']):
                    self._check_field(impact, 'assumption', str)
                    self._check_field(impact, 'impact', str)
                    self._check_field(impact, 'mitigation', str)
            else:
                self.errors.append("Missing 'assumptions.violations_impact' field")
        else:
            self.errors.append("Missing 'assumptions' field")

        # Implementation guidance
        if 'implementation_guidance' in handoff:
            guidance = handoff['implementation_guidance']
            self._check_field(guidance, 'library', str)
            self._check_field(guidance, 'function', str)
            self._check_field(guidance, 'parameters', dict)
            self._check_field(guidance, 'code_example', str)
        else:
            self.errors.append("Missing 'implementation_guidance' field")

        # Power analysis (optional)
        if 'power_analysis' in handoff:
            power = handoff['power_analysis']
            self._check_field(power, 'effect_size', (int, float))
            self._check_field(power, 'alpha', (int, float))
            self._check_field(power, 'power', (int, float))
            self._check_field(power, 'required_n', int)
            self._check_field(power, 'calculation_method', str)

        # Validation criteria
        if 'validation_criteria' in handoff:
            validation = handoff['validation_criteria']

            if 'diagnostic_checks' in validation:
                for i, check in enumerate(validation['diagnostic_checks']):
                    self._check_field(check, 'name', str)
                    self._check_field(check, 'method', str)
                    self._check_field(check, 'threshold', str)
            else:
                self.errors.append("Missing 'validation_criteria.diagnostic_checks' field")

            self._check_field(validation, 'sensitivity_analyses', list)
        else:
            self.errors.append("Missing 'validation_criteria' field")

        # MCMC config (optional)
        if 'mcmc_config' in handoff:
            mcmc = handoff['mcmc_config']
            self._check_field(mcmc, 'n_chains', int)
            self._check_field(mcmc, 'warmup', int)
            self._check_field(mcmc, 'samples', int)
            self._check_field(mcmc, 'thinning', int)

            if 'convergence_criteria' in mcmc:
                conv = mcmc['convergence_criteria']
                self._check_field(conv, 'ess_threshold', int)
                self._check_field(conv, 'rhat_threshold', (int, float))
            else:
                self.errors.append("Missing 'mcmc_config.convergence_criteria' field")

        # Interpretation guide
        if 'interpretation_guide' in handoff:
            interp = handoff['interpretation_guide']
            self._check_field(interp, 'result_format', str)
            self._check_field(interp, 'significant_threshold', (int, float))
            self._check_field(interp, 'interpretation_template', str)
        else:
            self.errors.append("Missing 'interpretation_guide' field")

        return self.errors


class CodeHandoffValidator(HandoffValidator):
    """Validator for developer -> code_review (code_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Task
        if 'task' in handoff:
            task = handoff['task']
            self._check_field(task, 'id', str)
            self._check_field(task, 'description', str)
            self._check_field(task, 'assigned_to', str)
        else:
            self.errors.append("Missing 'task' field")

        # Changes
        if 'changes' in handoff:
            changes = handoff['changes']

            if 'files_changed' in changes:
                for i, file in enumerate(changes['files_changed']):
                    self._check_field(file, 'path', str)
                    self._check_field(file, 'type', str)
                    if 'type' in file and file['type'] not in ['added', 'modified', 'deleted']:
                        self.errors.append(f"Invalid file change type: {file['type']}")
                    self._check_field(file, 'changes', str)
            else:
                self.errors.append("Missing 'changes.files_changed' field")

            self._check_field(changes, 'lines_added', int)
            self._check_field(changes, 'lines_removed', int)
        else:
            self.errors.append("Missing 'changes' field")

        # Summary (min 100 chars)
        if 'summary' in handoff:
            if not isinstance(handoff['summary'], str):
                self.errors.append(f"Invalid type for summary: expected str")
            elif len(handoff['summary']) < 100:
                self.errors.append(f"Summary too short: {len(handoff['summary'])} chars (min 100)")
        else:
            self.errors.append("Missing 'summary' field")

        # Test coverage
        if 'test_coverage' in handoff:
            coverage = handoff['test_coverage']
            self._check_field(coverage, 'new_lines', int)
            self._check_field(coverage, 'covered_lines', int)
            self._check_field(coverage, 'coverage_percent', (int, float))
        else:
            self.errors.append("Missing 'test_coverage' field")

        # Self-review checklist
        if 'self_review_checklist' in handoff:
            checklist = handoff['self_review_checklist']
            self._check_field(checklist, 'tests_pass', bool)
            self._check_field(checklist, 'ruff_clean', bool)
            self._check_field(checklist, 'mypy_clean', bool)
            self._check_field(checklist, 'documentation_updated', bool)
            self._check_field(checklist, 'type_hints_present', bool)
        else:
            self.errors.append("Missing 'self_review_checklist' field")

        self._check_field(handoff, 'open_questions', list)
        self._check_field(handoff, 'known_limitations', list)

        # Revision (optional)
        if 'revision' in handoff:
            revision = handoff['revision']
            self._check_field(revision, 'number', int)
            self._check_field(revision, 'changes_made', list)
            self._check_field(revision, 'previous_feedback_addressed', list)

        return self.errors


class ReviewHandoffValidator(HandoffValidator):
    """Validator for code_review -> merge (review_handoff)."""

    def validate(self) -> List[str]:
        super().validate()

        if 'handoff' not in self.data:
            return self.errors

        handoff = self.data['handoff']

        # Review
        if 'review' in handoff:
            review = handoff['review']
            self._check_field(review, 'reviewer', str)
            self._check_field(review, 'reviewed_date', str)
            if 'reviewed_date' in review:
                try:
                    datetime.fromisoformat(review['reviewed_date'].replace('Z', '+00:00'))
                except (ValueError, AttributeError):
                    self.errors.append(f"Invalid ISO8601 reviewed_date")

            self._check_field(review, 'status', str)
            if 'status' in review and review['status'] not in ['approved', 'changes_requested']:
                self.errors.append(f"Invalid review status: {review['status']}")
        else:
            self.errors.append("Missing 'review' field")

        # Automated checks
        if 'automated_checks' in handoff:
            checks = handoff['automated_checks']
            for check_name in ['ruff', 'mypy', 'tests']:
                self._check_field(checks, check_name, str)
                if check_name in checks and checks[check_name] not in ['pass', 'fail']:
                    self.errors.append(f"Invalid {check_name} status: {checks[check_name]}")
            self._check_field(checks, 'coverage', (int, float))
        else:
            self.errors.append("Missing 'automated_checks' field")

        # Manual review
        if 'manual_review' in handoff:
            manual = handoff['manual_review']
            for aspect in ['code_quality', 'documentation', 'testing', 'architecture']:
                self._check_field(manual, aspect, str)
                if aspect in manual and manual[aspect] not in ['pass', 'issues']:
                    self.errors.append(f"Invalid {aspect} status: {manual[aspect]}")
        else:
            self.errors.append("Missing 'manual_review' field")

        self._check_field(handoff, 'required_changes', list)
        self._check_field(handoff, 'suggestions', list)

        # Approval
        if 'approval' in handoff:
            approval = handoff['approval']
            self._check_field(approval, 'approved', bool)
            self._check_field(approval, 'conditions', list)
        else:
            self.errors.append("Missing 'approval' field")

        return self.errors


# Map handoff type to validator class
VALIDATORS = {
    'session_handoff': SessionHandoffValidator,
    'requirements_handoff': RequirementsHandoffValidator,
    'premortem_handoff': PremortемHandoffValidator,
    'architecture_handoff': ArchitectureHandoffValidator,
    'math_handoff': MathHandoffValidator,
    'stats_handoff': StatsHandoffValidator,
    'code_handoff': CodeHandoffValidator,
    'review_handoff': ReviewHandoffValidator,
}


def validate_handoff(handoff_path: Path, handoff_type: str) -> List[str]:
    """Validate handoff file against schema.

    Args:
        handoff_path: Path to YAML handoff file
        handoff_type: One of the 8 handoff types (session_handoff, requirements_handoff, etc.)

    Returns:
        List of error messages (empty if validation successful)
    """
    if handoff_type not in VALIDATORS:
        return [f"Unknown handoff type: {handoff_type}. Valid types: {', '.join(VALIDATORS.keys())}"]

    if not handoff_path.exists():
        return [f"Handoff file not found: {handoff_path}"]

    try:
        with open(handoff_path) as f:
            handoff_data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        return [f"Invalid YAML: {e}"]
    except Exception as e:
        return [f"Error reading file: {e}"]

    validator_class = VALIDATORS[handoff_type]
    validator = validator_class(handoff_data, handoff_path)
    return validator.validate()


def main():
    """CLI entry point."""
    if len(sys.argv) != 3:
        print("Usage: python3 validate-handoff.py <handoff_path> <handoff_type>", file=sys.stderr)
        print(f"\nValid handoff types: {', '.join(VALIDATORS.keys())}", file=sys.stderr)
        print("\nExample:", file=sys.stderr)
        print("  python3 validate-handoff.py /tmp/session/handoffs/phase1-requirements-handoff.yaml requirements_handoff", file=sys.stderr)
        sys.exit(2)

    handoff_path = Path(sys.argv[1])
    handoff_type = sys.argv[2]

    errors = validate_handoff(handoff_path, handoff_type)

    if errors:
        print(f"❌ Handoff validation FAILED: {handoff_path}", file=sys.stderr)
        print(f"\nErrors found ({len(errors)}):", file=sys.stderr)
        for i, error in enumerate(errors, 1):
            print(f"  {i}. {error}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"✅ Handoff validation PASSED: {handoff_path}")
        sys.exit(0)


if __name__ == '__main__':
    main()
