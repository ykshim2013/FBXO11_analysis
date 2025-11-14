# Rosetta Stability Predictions

> **⚠️ Data Files Status:** CSV data files to be added before publication. This README describes the expected data format and structure.

## Overview

This directory will contain Rosetta cartesian_ddg ΔΔG predictions for 45 FBXO11 missense variants across 270 structural contexts.

## Dataset Design

- **45 variants**: 24 pathogenic, 21 benign (wildtype excluded - see below)
- **6 contexts per variant**: 3 conformational seeds × 2 biological states (monomer, SKP1-bound)
- **10 technical replicates per context**: For statistical robustness
- **Total systems**: 270 (45 variants × 6 contexts)

## File Structure

```
rosetta/
├── aggregated_results.csv      # All 270 systems × 10 replicates
├── by_variant_summary.csv      # Per-variant statistics (mean ± SD)
└── README.md                   # This file
```

## Why Wildtype is Excluded

Rosetta cartesian_ddg calculates ΔΔG for mutations. WT→WT is mathematically undefined (no mutation), so wildtype systems are not included. This is **by design**, not a failure.

Result: 240/240 valid systems = **100% of analyzable variants**

See [`Documentation/WHY_ROSETTA_NOT_100_PERCENT.md`](../../Documentation/WHY_ROSETTA_NOT_100_PERCENT.md) for details.

## Data Columns

### `aggregated_results.csv`
- `variant`: Amino acid change (e.g., K135R)
- `classification`: pathogenic or benign
- `system_type`: monomer or complex
- `seed`: AlphaFold3 conformational seed (seed1, seed2, seed3)
- `replicate`: Technical replicate number (1-10)
- `ddg`: ΔΔG value in kcal/mol

### `by_variant_summary.csv`
- `variant`: Amino acid change
- `classification`: pathogenic or benign
- `ddg_mean`: Mean ΔΔG across all contexts
- `ddg_sd`: Standard deviation
- `ddg_median`: Median ΔΔG (for reference)
- `ddg_iqr`: Interquartile range (for reference)

## ΔΔG Interpretation

**Thresholds:**
- ΔΔG > 2.0 kcal/mol: Strongly destabilizing (likely pathogenic)
- 1.0 < ΔΔG ≤ 2.0: Moderately destabilizing (VUS/pathogenic)
- -1.0 ≤ ΔΔG ≤ 1.0: Neutral (VUS)
- ΔΔG < -1.0: Stabilizing (likely benign)

**Positive ΔΔG**: Mutation destabilizes protein (increases free energy)
**Negative ΔΔG**: Mutation stabilizes protein (decreases free energy)

## Key Results

| Classification | Mean ΔΔG | SD | p-value |
|----------------|---------|-------|---------|
| Pathogenic | 7.40 kcal/mol | 9.36 | 8.7×10⁻¹² |
| Benign | 1.31 kcal/mol | 3.30 | - |

**Welch's t-test:** p = 8.7×10⁻¹², Cohen's d = 0.84

## Citation

**Rosetta:**
Park H, Bradley P, Greisen P Jr, Liu Y, Mulligan VK, Kim DE, Baker D, DiMaio F (2016). Simultaneous Optimization of Biomolecular Energy Functions on Features from Small Molecules and Macromolecules. *J Chem Theory Comput*, 12(12):6201-6212.

Alford RF, Leaver-Fay A, Jeliazkov JR, et al. (2017). The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design. *J Chem Theory Comput*, 13(6):3031-3048.

## Usage

```python
import pandas as pd

# Load aggregated data
rosetta_all = pd.read_csv('rosetta/aggregated_results.csv')

# Load variant summaries
rosetta_summary = pd.read_csv('rosetta/by_variant_summary.csv')

# Filter pathogenic variants
pathogenic = rosetta_summary[rosetta_summary['classification'] == 'pathogenic']

# Calculate mean destabilization
print(f"Pathogenic mean: {pathogenic['ddg_mean'].mean():.2f} kcal/mol")
```

## Notes

- **Wildtype excluded by design**: WT→WT undefined for mutation-based ΔΔG
- **Parametric statistics**: Mean ± SD used (Rosetta output approximately normal)
- **Higher magnitude than FoldX**: Rosetta typically predicts larger ΔΔG values
- **Multi-conformational**: Captures conformational variability via AlphaFold3 ensemble

---

**Last updated**: 2025-11-14
