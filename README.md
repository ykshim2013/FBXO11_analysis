# FBXO11 Variant Stability Analysis

[![DOI](https://img.shields.io/badge/DOI-pending-blue)](https://github.com/ykshim2013/FBXO11_analysis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Computational framework for analyzing structural destabilization in *FBXO11* missense variants**

This repository contains computational analysis scripts, data, and protocols for the manuscript:

> **Structural Destabilization as a Pathogenic Mechanism in *FBXO11*-Related Neurodevelopmental Disorder: A Computational Framework for Variant Interpretation**
>
> Youngkyu Shim, M.D.
> Department of Pediatrics, Korea University Ansan Hospital

---

## Overview

*FBXO11*-related neurodevelopmental disorder is a rare Mendelian condition caused by pathogenic variants in the *FBXO11* gene, which encodes an F-box protein component of the SCF E3 ubiquitin ligase complex. This repository provides a comprehensive computational framework for analyzing the structural stability effects of *FBXO11* missense variants using:

- **FoldX 5.0** - Empirical force field-based stability predictions
- **Rosetta 2023.49** - Physics-based cartesian_ddg calculations
- **GROMACS 2024.4** - Molecular dynamics simulations
- **AlphaFold3** - Structural modeling

### Key Findings

- **45 variants analyzed**: 24 pathogenic, 21 benign (near-comprehensive coverage)
- **270 structural contexts**: Multi-conformational analysis across monomer and SKP1-bound states
- **300 ns MD simulations**: Revealing dynamic instability in pathogenic variants
- **Independence of stability and binding**: Pathogenic variants retain SKP1 binding despite destabilization

---

## Repository Structure

```
FBXO11_analysis/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── CITATION.cff                       # Citation information
│
├── data/                              # Computational data
│   ├── foldx/                        # FoldX ΔΔG predictions
│   │   ├── aggregated_results.csv    # All FoldX results (276 systems)
│   │   ├── by_variant_summary.csv    # Per-variant statistics
│   │   └── README.md                 # Data description
│   │
│   ├── rosetta/                      # Rosetta ΔΔG predictions
│   │   ├── aggregated_results.csv    # All Rosetta results (270 systems)
│   │   ├── by_variant_summary.csv    # Per-variant statistics
│   │   └── README.md                 # Data description
│   │
│   ├── md_trajectories/              # Molecular dynamics data
│   │   ├── trajectory_summary.csv    # RMSD, RMSF, SASA metrics
│   │   ├── md_analysis_results.csv   # Comprehensive MD analysis
│   │   ├── README.md                 # MD data description
│   │   └── download_instructions.md  # Large trajectory file access
│   │
│   └── alphafold3_structures/        # AlphaFold3 model information
│       ├── structure_metadata.csv    # pLDDT scores, RMSD values
│       └── README.md                 # Structure description
│
├── scripts/                          # Analysis scripts
│   ├── foldx_analysis/              # FoldX workflow
│   │   ├── 01_structure_preparation.py
│   │   ├── 02_foldx_batch_analysis.py
│   │   └── 03_aggregate_results.py
│   │
│   ├── rosetta_analysis/            # Rosetta workflow
│   │   ├── 01_prepare_structures.py
│   │   ├── 02_rosetta_ddg.py
│   │   └── 03_aggregate_results.py
│   │
│   ├── md_simulations/              # MD workflow
│   │   ├── 01_prepare_systems.py
│   │   ├── 02_run_md.sh
│   │   ├── 03_analyze_trajectories.py
│   │   └── 04_plot_results.py
│   │
│   ├── statistical_analysis/         # Statistical comparisons
│   │   ├── pathogenic_vs_benign.py
│   │   ├── roc_analysis.py
│   │   └── correlation_analysis.py
│   │
│   └── visualization/                # Figure generation
│       ├── figure1_stability_plots.py
│       ├── figure2_roc_curves.py
│       ├── figure3_md_rmsd.py
│       └── figure4_rmsf_profiles.py
│
├── analysis_outputs/                 # Generated results
│   ├── figures/                     # Publication figures
│   ├── tables/                      # Summary tables
│   └── statistics/                  # Statistical test results
│
└── docs/                            # Documentation
    ├── METHODS.md                   # Detailed computational protocols
    ├── DATA_DICTIONARY.md           # Data format descriptions
    ├── INSTALLATION.md              # Software installation guide
    └── TUTORIAL.md                  # Usage examples

```

---

## Requirements

### Software Dependencies

**Stability Prediction:**
- FoldX 5.0 (academic license: [FoldX website](https://foldxsuite.crg.eu/))
- Rosetta 2023.49 (academic license: [Rosetta Commons](https://www.rosettacommons.org/))

**Molecular Dynamics:**
- GROMACS 2024.4 or later ([GROMACS website](http://www.gromacs.org/))
- GPU support recommended (CUDA-enabled)

**Structure Prediction:**
- AlphaFold3 Server access ([AlphaFold Server](https://alphafoldserver.com))

**Python Environment:**
```bash
# Python 3.8+ required
python>=3.8
numpy>=1.24
pandas>=2.0
scipy>=1.11
matplotlib>=3.7
seaborn>=0.12
scikit-learn>=1.3
biopython>=1.81
MDAnalysis>=2.4
```

### Installation

```bash
# Clone repository
git clone https://github.com/ykshim2013/FBXO11_analysis.git
cd FBXO11_analysis

# Create conda environment
conda create -n fbxo11 python=3.10
conda activate fbxo11

# Install Python dependencies
pip install -r requirements.txt

# Install GROMACS (if not already installed)
# See docs/INSTALLATION.md for detailed instructions
```

---

## Usage

### Quick Start

```bash
# Activate environment
conda activate fbxo11

# Run statistical analysis on pre-computed data
cd scripts/statistical_analysis
python pathogenic_vs_benign.py

# Generate publication figures
cd ../visualization
python figure1_stability_plots.py
```

### Reproducing the Analysis

**1. FoldX Stability Predictions**
```bash
cd scripts/foldx_analysis
python 01_structure_preparation.py --input ../../data/alphafold3_structures/
python 02_foldx_batch_analysis.py --structures prepared_structures/
python 03_aggregate_results.py --output ../../data/foldx/
```

**2. Rosetta Stability Predictions**
```bash
cd scripts/rosetta_analysis
python 01_prepare_structures.py --input ../../data/alphafold3_structures/
python 02_rosetta_ddg.py --structures prepared_structures/
python 03_aggregate_results.py --output ../../data/rosetta/
```

**3. Molecular Dynamics Simulations**
```bash
cd scripts/md_simulations
python 01_prepare_systems.py --variants G421R,D910G,R138G,I300V,WT
bash 02_run_md.sh  # Requires GPU, ~48 hours for 15 simulations
python 03_analyze_trajectories.py --input ../../data/md_trajectories/
python 04_plot_results.py
```

**4. Statistical Analysis**
```bash
cd scripts/statistical_analysis
python pathogenic_vs_benign.py  # Mann-Whitney U, Welch's t-test
python roc_analysis.py          # ROC/PR curves, AUC calculation
python correlation_analysis.py  # FoldX-Rosetta correlation
```

See [`docs/TUTORIAL.md`](docs/TUTORIAL.md) for detailed usage examples.

---

## Data Availability

### Included in Repository

- **FoldX predictions**: Complete ΔΔG values for 276 systems (45 variants × 6 contexts)
- **Rosetta predictions**: Complete ΔΔG values for 270 systems
- **MD summary metrics**: RMSD, RMSF, SASA for 15 simulations (5 variants × 3 replicates)
- **AlphaFold3 metadata**: Structure quality metrics (pLDDT scores)

### Large Files (External Storage)

MD trajectory files (`.xtc`, `.tpr`, `.edr`) total ~55 GB and are available via:
- **Zenodo**: [DOI to be assigned upon publication]
- **Google Drive**: [Link in `data/md_trajectories/download_instructions.md`]

Each trajectory includes:
- `md.xtc`: Coordinates (3.6 GB, 2001 frames)
- `md.tpr`: Topology (16 MB)
- `md.edr`: Energy data (1.4 MB)
- `md.log`: Simulation log (1.3 MB)
- `md.gro`: Final structure (35 MB)

---

## Key Results

### Stability Discrimination

| Method | Pathogenic ΔΔG | Benign ΔΔG | p-value | Cohen's d | AUC |
|--------|---------------|-----------|---------|-----------|-----|
| **FoldX** | 1.88 ± 4.95 | 0.44 ± 2.79 | 1.0×10⁻⁵ | 0.49 | 0.65 |
| **Rosetta** | 7.40 ± 9.36 | 1.31 ± 3.30 | 8.7×10⁻¹² | 0.84 | 0.69 |
| **Combined** | - | - | - | - | **0.70** |

### Molecular Dynamics Results

| Variant | Type | RMSD (nm) | RMSF Fold Change | SKP1 Binding |
|---------|------|-----------|------------------|--------------|
| **G421R** | Pathogenic | 3.60 ± 2.04 | 1.10× | Normal |
| **D910G** | Pathogenic | 3.56 ± 0.55 | 1.20× | Normal |
| **R138G** | Pathogenic | 3.72 ± 1.61 | 1.05× | Normal |
| **I300V** | Benign | 2.86 ± 2.09 | 0.86× | Normal |
| **WT** | Reference | 2.61 ± 1.48 | 1.00× | Normal |

**Key finding**: D910G and R138G showed dynamic instability despite low static predictions, revealing limitations of FoldX/Rosetta for dynamic mechanisms.

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{shim2025fbxo11,
  title={Structural Destabilization as a Pathogenic Mechanism in FBXO11-Related Neurodevelopmental Disorder: A Computational Framework for Variant Interpretation},
  author={Shim, Youngkyu},
  journal={Frontiers in Genetics},
  year={2025},
  volume={TBD},
  pages={TBD},
  doi={TBD}
}
```

For software used in this analysis:

**FoldX**: Delgado et al. (2019). FoldX 5.0. *Bioinformatics* 35(20):4168-4169.
**Rosetta**: Alford et al. (2017). The Rosetta All-Atom Energy Function. *J Chem Theory Comput* 13(6):3031-3048.
**GROMACS**: Abraham et al. (2015). GROMACS. *SoftwareX* 1:19-25.
**AlphaFold3**: Abramson et al. (2024). Accurate structure prediction. *Nature* 630:493-500.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Note**: FoldX and Rosetta require separate academic licenses. This repository provides analysis scripts only.

---

## Contact

**Youngkyu Shim, M.D.**
Department of Pediatrics
Korea University Ansan Hospital
Korea University College of Medicine

📧 Email: ykshim2013@gmail.com
🔗 GitHub: [@ykshim2013](https://github.com/ykshim2013)

For questions about:
- **Computational methods**: Open an issue in this repository
- **Clinical interpretation**: Contact corresponding author
- **Data access**: See `data/*/README.md` files

---

## Acknowledgments

This work was supported by Korea University Ansan Hospital (grant O2515271).

**Computational resources**: MD simulations performed on GPU-accelerated workstation.

**Software acknowledgments**: FoldX, Rosetta, GROMACS, AlphaFold3, Python scientific stack.

---

## Version History

- **v1.0.0** (2025-01-XX): Initial release with manuscript publication
  - Complete FoldX and Rosetta predictions for 45 variants
  - MD simulations for 5 representative variants
  - Statistical analysis scripts
  - Publication figure generation code

---

## Contributing

This repository primarily serves as a data and code archive for the published manuscript. However, we welcome:

- Bug reports and fixes
- Documentation improvements
- Extension to additional variants or genes
- Integration with other variant interpretation tools

Please open an issue or pull request for any contributions.

---

**Last updated**: 2025-11-14
