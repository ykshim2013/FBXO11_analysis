# Molecular Dynamics Simulation Trajectories

> **⚠️ Data Files Status:** Summary CSV files to be added before publication. Full trajectory files (~55 GB) available via external storage (see download_instructions.md). This README describes the expected data format and structure.

## Overview

This directory will contain GROMACS 2024.4 molecular dynamics simulation summary data for 5 representative FBXO11 variants analyzed over 300 ns total simulation time.

## Dataset Design

- **5 variants**: 3 pathogenic (D910G, G421R, R138G), 1 benign (I300V), 1 wildtype (WT)
- **3 independent replicates per variant**: For statistical robustness
- **20 ns per replicate**: 60 ns total per variant
- **Total simulation time**: 300 ns (15 systems × 20 ns)

## File Structure

```
md_trajectories/
├── trajectory_summary.csv      # RMSD, RMSF, SASA metrics for all systems
├── md_analysis_results.csv     # Comprehensive MD analysis
├── README.md                   # This file
└── download_instructions.md    # How to access full trajectory files (~55 GB)
```

## Large Trajectory Files

**Note:** Full MD trajectory files (`.xtc`, `.tpr`, `.edr`) total ~55 GB and are **not included** in this GitHub repository. They are available via external storage:

- **Zenodo**: [DOI to be assigned upon publication]
- **Google Drive**: See [`download_instructions.md`](download_instructions.md)

Each trajectory replicate includes:
- `md.xtc`: Coordinates (3.6 GB, 2001 frames)
- `md.tpr`: Topology (16 MB)
- `md.edr`: Energy data (1.4 MB)
- `md.log`: Simulation log (1.3 MB)
- `md.gro`: Final structure (35 MB)

## Summary Data Columns

### `trajectory_summary.csv`
- `variant`: Amino acid change or WT
- `classification`: pathogenic, benign, or wildtype
- `replicate`: Independent replicate number (1-3)
- `rmsd_mean`: Mean backbone RMSD (nm)
- `rmsd_sd`: RMSD standard deviation
- `rmsd_max`: Maximum RMSD observed
- `rmsf_mean`: Mean per-residue flexibility (nm)
- `sasa_mean`: Mean solvent-accessible surface area (nm²)
- `rg_mean`: Mean radius of gyration (nm)

### `md_analysis_results.csv`
- Extended metrics including secondary structure analysis, hydrogen bonding, and energy decomposition

## Key Results

| Variant | Type | RMSD (nm) | RMSF Fold Change | Interpretation |
|---------|------|-----------|------------------|----------------|
| **G421R** | Pathogenic | 3.60 ± 2.04 | 1.10× | Extreme fluctuations (peak 7.2 nm) |
| **D910G** | Pathogenic | 3.56 ± 0.55 | 1.20× | Dynamic instability despite low static ΔΔG |
| **R138G** | Pathogenic | 3.72 ± 1.61 | 1.05× | Progressive drift (1.91→5.06 nm) |
| **I300V** | Benign | 2.86 ± 2.09 | 0.86× | Similar to wildtype |
| **WT** | Reference | 2.61 ± 1.48 | 1.00× | Baseline stability |

**Key finding:** D910G and R138G showed dynamic instability despite FoldX/Rosetta predictions below pathogenic thresholds, revealing fundamental limitations of static stability methods.

## Simulation Protocol

**Software:** GROMACS 2024.4 with GPU acceleration
**Force field:** AMBER99SB-ILDN
**Water model:** TIP3P
**Ensemble:** NPT (1 bar, 310 K)
**Integration:** 2 fs timestep
**Constraints:** LINCS algorithm for bonds
**Electrostatics:** PME with 1.0 nm cutoff
**Output frequency:** 10 ps (2001 frames per 20 ns)

Full protocol: [`../../docs/METHODS.md`](../../docs/METHODS.md)

## Citation

**GROMACS:**
Abraham MJ, Murtola T, Schulz R, Páll S, Smith JC, Hess B, Lindahl E (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX*, 1:19-25.

**AMBER Force Field:**
Lindorff-Larsen K, Piana S, Palmo K, Maragakis P, Klepeis JL, Dror RO, Shaw DE (2010). Improved side-chain torsion potentials for the Amber ff99SB protein force field. *Proteins*, 78(8):1950-8.

## Usage

```python
import pandas as pd
import numpy as np

# Load summary data
md_summary = pd.read_csv('md_trajectories/trajectory_summary.csv')

# Calculate mean RMSD for pathogenic variants
pathogenic = md_summary[md_summary['classification'] == 'pathogenic']
print(f"Pathogenic mean RMSD: {pathogenic['rmsd_mean'].mean():.2f} nm")

# Load full trajectory (if downloaded)
# import MDAnalysis as mda
# u = mda.Universe('trajectories/D910G_rep1/md.tpr', 'trajectories/D910G_rep1/md.xtc')
```

## Limitations

**Important:** Only 5 of 45 variants (11%) have MD simulation data. Comprehensive MD profiling across all variants is needed to determine the true prevalence of dynamic instability mechanisms.

D910G and R138G cases demonstrate that static methods (FoldX/Rosetta) may miss pathogenic variants operating through dynamic conformational mechanisms rather than thermodynamic destabilization.

## Future Work

Extended simulations (100-500 ns or microsecond-scale enhanced sampling) recommended to:
1. Capture slow-onset conformational transitions
2. Validate static predictions for all 45 variants
3. Identify additional dynamic instability mechanisms

---

**Last updated**: 2025-11-14
