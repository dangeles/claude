# Rasmussen's Abstraction Hierarchy

Framework for systematic problem abstraction to enable cross-domain analogical reasoning.

## Overview

Rasmussen's Abstraction Hierarchy (1985) provides a five-level framework for describing systems at different levels of abstraction. Originally developed for process control and human-system interaction, it's highly effective for cross-domain knowledge transfer.

## The Five Levels

```
Level 1: Functional Purpose (Why?)
    ↓
Level 2: Abstract Function (What flows?)
    ↓
Level 3: Generalized Function (What type?)
    ↓
Level 4: Physical Function (How?)
    ↓
Level 5: Physical Form (What?)
```

### Level 1: Functional Purpose

**Question**: Why does the system exist? What goal does it serve?

**Characteristics**:
- Highest abstraction
- Domain-independent
- Answers "what value does this create?"

**Example for B2B SaaS Churn**:
- "Maintain ongoing value exchange relationships"
- "Sustain voluntary associations over time"

### Level 2: Abstract Function

**Question**: What flows through the system? What is being transformed?

**Characteristics**:
- Identifies flows (information, materials, energy)
- Defines transformations
- Maps constraints

**Example for B2B SaaS Churn**:
- "Retention through continuous value delivery"
- Flow: Value perception → Decision → Action
- Constraint: Customer has exit right

### Level 3: Generalized Function

**Question**: What category of process is this?

**Characteristics**:
- Process categorization
- Standard patterns
- Discipline-independent

**Example for B2B SaaS Churn**:
- "Relationship maintenance under optional commitment"
- Category: Voluntary retention process

### Level 4: Physical Function

**Question**: How does the specific mechanism work?

**Characteristics**:
- Domain-specific
- Concrete mechanisms
- Implementation details matter

**Example for B2B SaaS Churn**:
- "Subscription renewal decision-making"
- Mechanism: Monthly evaluation → Cancellation option → Renewal

### Level 5: Physical Form

**Question**: What does the concrete implementation look like?

**Characteristics**:
- Most concrete
- Specific to domain
- Observable features

**Example for B2B SaaS Churn**:
- "Monthly SaaS billing with one-click cancellation"
- Form: Stripe subscriptions, email notifications

## Using the Hierarchy for Analogies

### Ascend to Find Commonality

To find analogous domains, ascend the hierarchy:

```
Physical Form (SaaS billing)
    ↑
Physical Function (Subscription decisions)
    ↑
Generalized Function (Voluntary retention) ← Look for analogies here
    ↑
Abstract Function (Value exchange)
    ↑
Functional Purpose (Maintain relationships)
```

At **Generalized Function** or **Abstract Function**, ask:
"What other domains have voluntary retention processes?"

Answer: Healthcare (patient retention), Fitness (gym membership), Professional services (client retention)

### Descend to Apply Solutions

Once analogous domain found, descend to implement:

```
Healthcare: "Continuity of Care Model"
    ↓
Generalized Function: Consistent relationship maintenance
    ↓
Physical Function: Provider-patient continuity
    ↓
Physical Form (Healthcare): Same doctor for each visit
    ↓
[TRANSFER]
    ↓
Physical Form (SaaS): Same CSM for each account
```

## Domain Distance and Abstraction Level

### Near Field (Low Abstraction)

Connect at Level 3-4:
- Example: SaaS → Consumer subscriptions
- Shared: Subscription model, voluntary churn
- Different: B2B vs B2C decision process

### Mid Field (Medium Abstraction)

Connect at Level 2-3:
- Example: SaaS → Healthcare patient retention
- Shared: Voluntary retention, value delivery
- Different: Regulatory context, switching costs

### Far Field (High Abstraction)

Connect at Level 1-2:
- Example: SaaS → Wildlife habitat conservation
- Shared: Voluntary association, retention through value
- Different: Commercial vs biological, human vs animal

**Rule**: The more distant the domains, the higher the abstraction needed to connect them.

## Example 1: Software Problem

### Original Problem
"Our machine learning model training is too slow"

### Rasmussen Hierarchy

| Level | Content |
|-------|---------|
| **1. Purpose** | Minimize time to deploy validated models |
| **2. Abstract Function** | Transformation: Data → Trained model. Constraint: Computational resources finite |
| **3. Generalized Function** | Optimization under resource constraints |
| **4. Physical Function** | Iterative gradient descent on GPU clusters |
| **5. Physical Form** | PyTorch on AWS EC2 p3.8xlarge instances |

### Analogous Domains (Level 3)

"Optimization under resource constraints" appears in:
- **Manufacturing** (production scheduling)
- **Logistics** (route optimization)
- **Finance** (portfolio optimization)
- **Biology** (metabolic efficiency)

### Transfer Example

From **Manufacturing (Near Field)**:
- Technique: Batch processing with work-in-progress limits
- Transfer: Mini-batch gradient descent with memory-aware batch sizing
- Adaptation: Replace "production units" with "training examples"

## Example 2: Hardware Problem

### Original Problem
"Our circuit overheats during peak load"

### Rasmussen Hierarchy

| Level | Content |
|-------|---------|
| **1. Purpose** | Maintain operational stability under variable demand |
| **2. Abstract Function** | Heat dissipation: Energy in → Heat out. Constraint: Thermal capacity limited |
| **3. Generalized Function** | Resource management under capacity limits |
| **4. Physical Function** | Conduction/convection heat transfer via heatsink |
| **5. Physical Form** | Aluminum heatsink with fan on CPU |

### Analogous Domains (Level 2-3)

"Resource management under capacity limits" and "Heat dissipation" appears in:
- **HVAC** (building cooling)
- **Automotive** (engine cooling)
- **Human physiology** (thermoregulation)
- **Data centers** (server cooling)

### Transfer Example

From **Human Physiology (Far Field)**:
- Technique: Sweating (evaporative cooling triggered by temperature sensors)
- Transfer: Phase-change cooling activated by temperature threshold
- Adaptation: Replace "sweat glands" with "heat pipes", "evaporation" with "phase-change material"
- Warning: Far field - requires validation of structural mapping

## Using Hierarchy in POV-Expansion

### Stage 2: Problem Abstraction

The pov-abstractor-classifier agent:
1. Takes user's problem (Level 5)
2. Ascends through levels asking key questions
3. Populates all 5 levels
4. Creates abstract representation (Levels 1-3 most important for analogies)

### Stage 3: Domain Classification

Domain distance = abstraction level needed:
- **Near**: Can connect at Level 3-4
- **Mid**: Need Level 2-3
- **Far**: Require Level 1-2

### Stage 4: Perspective Generation

Perspective analysts receive Level 1-3 to find structural analogies in their domains.

## Common Pitfalls

### Pitfall 1: Staying Too Concrete

❌ "Our SaaS has high churn" → "Look at other SaaS companies"
✅ "Our SaaS has high churn" → Abstract to "voluntary retention" → Find all domains with voluntary retention

### Pitfall 2: Abstracting Too Far

❌ "Reduce churn" → "Achieve goals" → Too vague to find analogies
✅ "Reduce churn" → "Maintain voluntary relationships" → Specific enough for structural mapping

### Pitfall 3: Skipping Levels

❌ Jump from Level 5 (SaaS billing) directly to Level 1 (relationships)
✅ Progress through levels systematically to maintain structural coherence

### Pitfall 4: Forced Analogies

❌ "SaaS is like a gym membership" (surface similarity: both subscriptions)
✅ Check: Do both have same constraints? Same value delivery mechanisms? Same churn dynamics?

## Quality Criteria

A good abstraction:
- [ ] All 5 levels populated with non-placeholder content
- [ ] Each level answers its key question
- [ ] Levels are coherent (each flows logically from previous)
- [ ] Domain-specific jargon removed at Levels 1-3
- [ ] Constraints and flows explicitly identified at Level 2
- [ ] Physical details preserved at Levels 4-5 for context

## Advanced: Multiple Abstractions

Complex problems may benefit from multiple abstraction paths:

**Problem**: "Reduce customer churn in SaaS"

**Abstraction Path A** (Relationship focus):
1. Purpose: Maintain voluntary relationships
2. Abstract Function: Value perception management
3. Generalized Function: Retention process

**Abstraction Path B** (Economic focus):
1. Purpose: Maximize customer lifetime value
2. Abstract Function: Revenue stream optimization
3. Generalized Function: Investment return process

Different paths lead to different analogous domains:
- Path A → Healthcare, education, professional services
- Path B → Finance, portfolio management, asset management

## Summary

Rasmussen's Abstraction Hierarchy enables:
- **Systematic abstraction**: From concrete to general
- **Cross-domain search**: Find analogies at appropriate abstraction level
- **Solution transfer**: Descend from abstract solution to concrete implementation
- **Distance assessment**: Measure how far domains are apart

**Key principle**: Higher abstraction = more distant but potentially more novel analogies.
