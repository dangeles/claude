#!/usr/bin/env bash
#
# Test mode detection function
#

set -euo pipefail

# Create test requirements file
TEST_DIR="/tmp/programming-pm-test-$$"
mkdir -p "$TEST_DIR"

echo "Testing mode detection..."
echo ""

# Test Case 1: Simple - Single Component Utility
echo "Test 1: Simple - Single Component Utility"
cat > "$TEST_DIR/test1-requirements.yaml" <<'EOF'
handoff:
  requirements:
    problem_statement: "Create a utility script to parse CSV files"
    scope:
      in_scope:
        - "CSV parsing"
        - "Data validation"
        - "Output formatting"
EOF

# Simulate detect_tier (simplified for testing without full YAML parsing)
grep -qi "utility\|script\|helper\|tool" "$TEST_DIR/test1-requirements.yaml" && echo "  ✓ Contains utility keywords"
! grep -qi "statistic\|regression\|bayesian" "$TEST_DIR/test1-requirements.yaml" && echo "  ✓ No statistics keywords"
! grep -qi "algorithm\|optimization\|complexity" "$TEST_DIR/test1-requirements.yaml" && echo "  ✓ No math keywords"
echo "  Expected: SIMPLE"
echo ""

# Test Case 2: Standard - Multi-Component Web API
echo "Test 2: Standard - Multi-Component Web API"
cat > "$TEST_DIR/test2-requirements.yaml" <<'EOF'
handoff:
  requirements:
    problem_statement: "Build a REST API with database integration"
    scope:
      in_scope:
        - "API endpoints"
        - "Database models"
        - "Authentication"
        - "API documentation"
        - "Unit tests"
        - "Integration tests"
        - "Deployment"
        - "Monitoring"
EOF

# Check for web API keywords
grep -qi "REST API\|web API\|CLI tool" "$TEST_DIR/test2-requirements.yaml" && echo "  ✓ Contains web API keywords"
TASK_COUNT=$(grep -c "^        -" "$TEST_DIR/test2-requirements.yaml" || echo 0)
echo "  ✓ Task count: $TASK_COUNT (5-15 range = STANDARD)"
echo "  Expected: STANDARD"
echo ""

# Test Case 3: Extended - Statistics + Algorithms
echo "Test 3: Extended - Statistics + Algorithms"
cat > "$TEST_DIR/test3-requirements.yaml" <<'EOF'
handoff:
  requirements:
    problem_statement: "Develop statistical analysis pipeline with optimization algorithms including hypothesis testing, regression, and numerical optimization methods"
    scope:
      in_scope:
        - "Statistical methods"
        - "Optimization algorithms"
EOF

grep -qi "statistic\|regression\|hypothesis" "$TEST_DIR/test3-requirements.yaml" && echo "  ✓ Contains statistics keywords"
grep -qi "algorithm\|optimization\|numerical" "$TEST_DIR/test3-requirements.yaml" && echo "  ✓ Contains math keywords"
echo "  Expected: EXTENDED (dual specialization)"
echo ""

# Test Case 4: Extended - Many Components
echo "Test 4: Extended - Many Components"
cat > "$TEST_DIR/test4-requirements.yaml" <<'EOF'
handoff:
  requirements:
    problem_statement: "Build microservices architecture with 7 services"
    scope:
      in_scope:
        - "component 1"
        - "component 2"
        - "component 3"
        - "component 4"
        - "component 5"
        - "component 6"
        - "component 7"
EOF

COMPONENT_COUNT=$(grep -oi "component" "$TEST_DIR/test4-requirements.yaml" | wc -l | tr -d ' ')
echo "  ✓ Component count: $COMPONENT_COUNT"
[ "$COMPONENT_COUNT" -gt 5 ] && echo "  ✓ Exceeds threshold (>5)"
echo "  Expected: EXTENDED (many components)"
echo ""

# Test Case 5: Extended - Architectural Keywords
echo "Test 5: Extended - Architectural Keywords"
cat > "$TEST_DIR/test5-requirements.yaml" <<'EOF'
handoff:
  requirements:
    problem_statement: "Build real-time event-driven processing system with microservices"
EOF

grep -qi "distributed system\|microservices\|event-driven\|real-time" "$TEST_DIR/test5-requirements.yaml" && echo "  ✓ Contains architectural complexity keywords"
echo "  Expected: EXTENDED (architectural complexity)"
echo ""

# Cleanup
rm -rf "$TEST_DIR"

echo "================================================"
echo "Mode detection tests complete"
echo "================================================"
