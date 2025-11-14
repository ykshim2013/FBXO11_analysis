# FoldX Stability Predictions

## Overview

This directory contains FoldX 5.0 ΔΔG predictions for 45 FBXO11 missense variants across 276 structural contexts.

## Dataset Design

- **45 variants**: 24 pathogenic, 21 benign
- **6 contexts per variant**: 3 conformational seeds × 2 biological states (monomer, SKP1-bound)
- **10 technical replicates per context**: For statistical robustness
- **Total systems**: 276 (45 variants × 6 contexts + wildtype)

## File Structure

```
foldx/
├── aggregated_results.csv      # All 276 systems × 10 replicates
├── by_variant_summary.csv      # Per-variant statistics (median ± IQR)
└── README.md                   # This file
```

## Data Columns

### `aggregated_results.csv`
- `variant`: Amino acid change (e.g., K135R)
- `classification`: pathogenic, benign, or wildtype
- `system_type`: monomer or complex
- `seed`: AlphaFold3 conformational seed (seed1, seed2, seed3)
- `replicate`: Technical replicate number (1-10)
- `ddg`: ΔΔG value in kcal/mol

### `by_variant_summary.csv`
- `variant`: Amino acid change
- `classification`: pathogenic, benign, or wildtype
- `ddg_median`: Median ΔΔG across all contexts
- `ddg_iqr`: Interquartile range
- `ddg_mean`: Mean ΔΔG (for reference)
- `ddg_sd`: Standard deviation (for reference)

## ΔΔG Interpretation

**Thresholds:**
- ΔΔG > 1.0 kcal/mol: Destabilizing (likely pathogenic)
- -1.0 ≤ ΔΔG ≤ 1.0: Neutral (VUS)
- ΔΔG < -1.0: Stabilizing (likely benign)

**Positive ΔΔG**: Mutation destabilizes protein (increases free energy)
**Negative ΔΔG**: Mutation stabilizes protein (decreases free energy)

## Key Results

| Classification | Median ΔΔG | IQR | p-value |
|----------------|-----------|-----|---------|
| Pathogenic | 1.88 kcal/mol | 4.95 | 1.0×10⁻⁵ |
| Benign | 0.44 kcal/mol | 2.79 | - |

**Mann-Whitney U test:** p = 1.0×10⁻⁵, Cohen's d = 0.49

## Citation

**FoldX:**
Delgado J, Radusky LG, Cianferoni D, Serrano L (2019). FoldX 5.0: working with RNA, small molecules and a new graphical interface. *Bioinformatics*, 35(20):4168-4169.

## Usage

```python
import pandas as pd

# Load aggregated data
foldx_all = pd.read_csv('foldx/aggregated_results.csv')

# Load variant summaries
foldx_summary = pd.read_csv('foldx/by_variant_summary.csv')

# Filter pathogenic variants
pathogenic = foldx_summary[foldx_summary['classification'] == 'pathogenic']

# Calculate mean destabilization
print(f"Pathogenic median: {pathogenic['ddg_median'].median():.2f} kcal/mol")
```

## Notes

- **Wildtype included**: WT analyzed across all 6 contexts for reference
- **Robust statistics**: Median ± IQR used due to outlier presence
- **Multi-conformational**: Captures conformational variability via AlphaFold3 ensemble

---

**Last updated**: 2025-11-14
