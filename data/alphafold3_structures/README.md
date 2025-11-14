# AlphaFold3 Structure Information

> **⚠️ Data Files Status:** Metadata CSV file to be added before publication. This README describes the expected data format and structure.

## Overview

This directory will contain metadata for AlphaFold3 structural models used in FoldX and Rosetta stability predictions.

## Structure Design

### Monomer Structures
- **Chain:** FBXO11 full-length (943 residues)
- **Seeds:** 3 conformational variants
  - seed1: rank_0 (highest pLDDT)
  - seed2: rank_3 (structural diversity)
  - seed3: rank_4 (structural diversity)
- **Structural diversity:** Pairwise RMSD > 19 Å between seeds

### Complex Structures (FBXO11-SKP1)
- **Chain A:** FBXO11 (943 residues)
- **Chain B:** SKP1 (163 residues)
- **Seeds:** 3 conformational variants
  - seed1: rank_001 (highest interface confidence)
  - seed2: rank_003 (structural diversity)
  - seed3: rank_005 (structural diversity)
- **Interface quality:** pLDDT > 80 for F-box domain-SKP1 interface

## File Structure

```
alphafold3_structures/
├── structure_metadata.csv    # pLDDT scores, RMSD values, quality metrics
└── README.md                 # This file
```

## Structure Metadata

The `structure_metadata.csv` file contains:
- `structure_id`: Unique identifier (e.g., monomer_seed1, complex_seed2)
- `pdb_file`: Original AlphaFold3 output filename
- `mean_plddt`: Average pLDDT score across all residues
- `f_box_plddt`: pLDDT score for F-box domain (residues 358-423)
- `interface_plddt`: pLDDT for FBXO11-SKP1 interface (complex only)
- `rmsd_to_seed1`: Backbone RMSD relative to seed1 (Å)

## Quality Metrics

### pLDDT Interpretation
- **> 90:** Very high confidence (typically correct)
- **70-90:** High confidence (generally reliable)
- **50-70:** Low confidence (may be unreliable)
- **< 50:** Very low confidence (likely incorrect)

**All seeds selected have pLDDT > 70 for F-box domain** to ensure reliable stability predictions.

### Structural Diversity
Multi-conformational analysis with RMSD > 19 Å between seeds captures:
- Conformational heterogeneity
- Domain flexibility
- Binding-induced structural changes

## AlphaFold3 Server Details

**Platform:** https://alphafoldserver.com
**Version:** AlphaFold3 (May 2024 release)
**Prediction date:** 2024-12-XX
**Models per input:** 5 ranked predictions

## Usage

Structure files used for:
1. FoldX stability predictions (RepairPDB + BuildModel)
2. Rosetta stability predictions (Relax + cartesian_ddg)
3. MD simulation starting structures (selected variants only)

## Citation

**AlphaFold3:**
Abramson J, Adler J, Dunger J, et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*, 630:493-500.

---

**Last updated**: 2025-11-14
