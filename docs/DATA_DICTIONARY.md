# Data Dictionary

Complete description of all data file formats and column definitions.

## Table of Contents

1. [FoldX Results](#foldx-results)
2. [Rosetta Results](#rosetta-results)
3. [MD Trajectory Data](#md-trajectory-data)
4. [AlphaFold3 Metadata](#alphafold3-metadata)

---

## FoldX Results

### `data/foldx/aggregated_results.csv`

Complete dataset with all systems and technical replicates.

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `variant` | string | Amino acid substitution | K135R |
| `classification` | string | Clinical classification | pathogenic, benign, wildtype |
| `system_type` | string | Biological context | monomer, complex |
| `seed` | string | AlphaFold3 conformational seed | seed1, seed2, seed3 |
| `replicate` | integer | Technical replicate number (1-10) | 5 |
| `ddg` | float | FoldX ΔΔG (kcal/mol) | 2.45 |
| `total_energy` | float | Total system energy (kcal/mol) | -823.12 |
| `backbone_hbond` | float | Backbone H-bond energy | -245.67 |
| `sidechain_hbond` | float | Sidechain H-bond energy | -89.34 |
| `van_der_waals` | float | Van der Waals energy | -412.89 |
| `electrostatics` | float | Electrostatic energy | -75.22 |

**Rows:** 2,760 (276 systems × 10 replicates)

**Sign convention:**
- Positive ΔΔG: Destabilizing mutation
- Negative ΔΔG: Stabilizing mutation

### `data/foldx/by_variant_summary.csv`

Per-variant statistics aggregated across all contexts.

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `variant` | string | Amino acid substitution | K135R |
| `classification` | string | Clinical classification | pathogenic |
| `ddg_median` | float | Median ΔΔG across all contexts | 3.12 |
| `ddg_iqr` | float | Interquartile range | 2.45 |
| `ddg_mean` | float | Mean ΔΔG (for reference) | 3.34 |
| `ddg_sd` | float | Standard deviation | 1.89 |
| `ddg_min` | float | Minimum ΔΔG | -0.23 |
| `ddg_max` | float | Maximum ΔΔG | 8.91 |
| `n_systems` | integer | Number of systems analyzed | 60 |

**Rows:** 45 (24 pathogenic + 21 benign + 1 wildtype)

**Note:** Median ± IQR used due to presence of outliers in FoldX predictions.

---

## Rosetta Results

### `data/rosetta/aggregated_results.csv`

Complete dataset with all systems and technical replicates.

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `variant` | string | Amino acid substitution | K135R |
| `classification` | string | Clinical classification | pathogenic, benign |
| `system_type` | string | Biological context | monomer, complex |
| `seed` | string | AlphaFold3 conformational seed | seed1, seed2, seed3 |
| `replicate` | integer | Technical replicate number (1-10) | 5 |
| `ddg` | float | Rosetta ΔΔG (kcal/mol) | 8.23 |
| `total_score` | float | Total Rosetta energy (REU) | -1245.67 |
| `fa_atr` | float | Lennard-Jones attractive | -823.45 |
| `fa_rep` | float | Lennard-Jones repulsive | 234.12 |
| `fa_sol` | float | Solvation energy | 156.78 |
| `fa_elec` | float | Coulombic electrostatics | -89.34 |
| `hbond_sc` | float | Sidechain H-bonds | -45.67 |

**Rows:** 2,700 (270 systems × 10 replicates)

**Note:** Wildtype excluded (WT→WT undefined for mutation-based ΔΔG)

**Sign convention:**
- Positive ΔΔG: Destabilizing mutation
- Negative ΔΔG: Stabilizing mutation

### `data/rosetta/by_variant_summary.csv`

Per-variant statistics aggregated across all contexts.

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `variant` | string | Amino acid substitution | K135R |
| `classification` | string | Clinical classification | pathogenic |
| `ddg_mean` | float | Mean ΔΔG across all contexts | 9.45 |
| `ddg_sd` | float | Standard deviation | 5.67 |
| `ddg_median` | float | Median ΔΔG (for reference) | 7.89 |
| `ddg_iqr` | float | Interquartile range | 6.23 |
| `ddg_min` | float | Minimum ΔΔG | 1.23 |
| `ddg_max` | float | Maximum ΔΔG | 24.56 |
| `n_systems` | integer | Number of systems analyzed | 60 |

**Rows:** 44 (24 pathogenic + 21 benign, wildtype excluded)

**Note:** Mean ± SD used (Rosetta predictions approximately normal).

---

## MD Trajectory Data

### `data/md_trajectories/trajectory_summary.csv`

Summary metrics for all MD simulations.

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `variant` | string | Amino acid substitution or WT | D910G |
| `classification` | string | Clinical classification | pathogenic |
| `replicate` | integer | Independent replicate (1-3) | 2 |
| `rmsd_mean` | float | Mean backbone RMSD (nm) | 3.56 |
| `rmsd_sd` | float | RMSD standard deviation (nm) | 0.55 |
| `rmsd_max` | float | Maximum RMSD (nm) | 4.89 |
| `rmsf_mean` | float | Mean per-residue flexibility (nm) | 0.18 |
| `sasa_mean` | float | Mean SASA (nm²) | 285.34 |
| `rg_mean` | float | Mean radius of gyration (nm) | 2.45 |
| `simulation_time` | integer | Total simulation time (ns) | 20 |

**Rows:** 15 (5 variants × 3 replicates)

**RMSD interpretation:**
- < 2.0 nm: Stable structure
- 2.0-4.0 nm: Moderate flexibility
- > 4.0 nm: High instability or unfolding

### `data/md_trajectories/md_analysis_results.csv`

Extended MD analysis with additional metrics.

| Column | Type | Description |
|--------|------|-------------|
| `variant` | string | Amino acid substitution |
| `classification` | string | Clinical classification |
| `replicate` | integer | Replicate number |
| `rmsd_mean` | float | Mean backbone RMSD (nm) |
| `rmsf_mean` | float | Mean per-residue RMSF (nm) |
| `sasa_mean` | float | Mean SASA (nm²) |
| `rg_mean` | float | Mean Rg (nm) |
| `helix_content` | float | α-helix percentage |
| `sheet_content` | float | β-sheet percentage |
| `hbond_intra_mean` | float | Intramolecular H-bonds |
| `hbond_water_mean` | float | Protein-water H-bonds |
| `energy_potential_mean` | float | Mean potential energy (kJ/mol) |
| `energy_kinetic_mean` | float | Mean kinetic energy (kJ/mol) |
| `temperature_mean` | float | Mean temperature (K) |
| `pressure_mean` | float | Mean pressure (bar) |

**Rows:** 15 (5 variants × 3 replicates)

---

## AlphaFold3 Metadata

### `data/alphafold3_structures/structure_metadata.csv`

Quality metrics for AlphaFold3 structural models.

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `structure_id` | string | Unique identifier | monomer_seed1 |
| `type` | string | Structure type | monomer, complex |
| `seed` | string | Seed identifier | seed1 |
| `pdb_file` | string | Original AlphaFold3 filename | fold_rank_0.pdb |
| `mean_plddt` | float | Average pLDDT score (0-100) | 78.5 |
| `f_box_plddt` | float | F-box domain pLDDT (residues 358-423) | 85.3 |
| `interface_plddt` | float | Interface pLDDT (complex only) | 82.1 |
| `rmsd_to_seed1` | float | Backbone RMSD to seed1 (Å) | 19.4 |
| `num_residues` | integer | Total residues | 943 |
| `num_models` | integer | AlphaFold3 models generated | 5 |
| `rank` | integer | Original AlphaFold3 rank | 0 |

**Rows:** 6 (3 monomer + 3 complex seeds)

**pLDDT interpretation:**
- > 90: Very high confidence
- 70-90: High confidence
- 50-70: Low confidence
- < 50: Very low confidence

---

## Data Types Reference

### String Types
- `variant`: Single letter amino acid code + position + new amino acid (e.g., K135R)
- `classification`: {pathogenic, benign, wildtype}
- `system_type`: {monomer, complex}
- `seed`: {seed1, seed2, seed3}

### Numeric Types
- `float`: Floating-point numbers (ΔΔG, energies, RMSD)
- `integer`: Whole numbers (replicate numbers, counts)

### Units
- **Energy:** kcal/mol (FoldX), REU (Rosetta Energy Units), kJ/mol (GROMACS)
- **Distance:** nm (nanometers) for MD, Å (Ångstroms) for structures
- **Time:** ns (nanoseconds)
- **Temperature:** K (Kelvin)
- **Pressure:** bar

---

## Missing Data Codes

- **NaN:** Not available or not applicable
- **Empty cell:** Data not collected for this system
- **Wildtype excluded (Rosetta):** WT→WT undefined for mutation-based ΔΔG

---

**Last updated**: 2025-11-14
