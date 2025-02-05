# Clause Weighting Parameters Explained

Definition: tuned weight = original clause weight / average clause weight

- Non-Partial
  - MaxSAT
    - Soft Clause Weight Initialization: `sWeightBoundCoeff` (397).
    - Soft Clause Weight Update:
      - Probabilistically smooth weights with probability `sSmoothProb` (200000 = 0.2%)
        by reducing weights of satisfied clauses by `sWeightInc` (1).
      - Otherwise, increase weights of falsified clauses by `sWeightInc` (1),
        bounded by `sWeightBoundCoeff` (397) + `sWeightBoundOffset` (550).
  - Weighted MaxSAT
    - Soft Clause Weight Initialization: `sWeightBoundCoeff` (1000) * tuned weight.
    - Soft Clause Weight Update:
      - Probabilistically smooth weights with probability `sSmoothProb` (100000 = 0.1%)
        by reducing weights of satisfied clauses by `sWeightInc` (3).
      - Otherwise, increase weights of falsified clauses by `sWeightInc` (3),
        bounded by `sWeightBoundCoeff` (1000) * tuned weight + `sWeightBoundOffset` (10).
- Partial
  - Partial MaxSAT
    - Hard Clause Weight Initialization: `hWeightInit` (1).
    - Soft Clause Weight Initialization: 0.
    - Hard Clause Weight Update: Increment weights of falsified clauses by `hWeightInc` (1).
    - Soft Clause Weight Update: If there are no falsified hard clauses,
      increase weights of all soft clauses by `sWeightInc` (1).
  - Weighted Partial MaxSAT
    - Hard Clause Weight Initialization: `hWeightInit` (1).
    - Soft Clause Weight Initialization: 0.
    - Hard Clause Weight Update: Increment weights of falsified clauses by `hWeightInc` (5).
    - Soft Clause Weight Update: If there are no falsified hard clauses,
      increase weights of all soft clauses by `sWeightInc` (1) * tuned weight.
