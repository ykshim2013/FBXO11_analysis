# MD Trajectory Download Instructions

## Overview

Full molecular dynamics trajectory files total approximately **55 GB** and are not included in this GitHub repository. They are available via external storage.

## Storage Locations

### Zenodo (Recommended)
**DOI:** [To be assigned upon publication]
**URL:** [https://zenodo.org/record/XXXXXX](https://zenodo.org/record/XXXXXX)

**Advantages:**
- Permanent DOI for citation
- Version control
- Academic standard for data sharing

### Google Drive (Alternative)
**URL:** [To be provided]

**Advantages:**
- Faster download speeds
- No download limits
- Direct file browsing

## File Organization

Each variant has 3 replicate directories with identical structure:

```
trajectories/
├── pathogenic/
│   ├── D910G/
│   │   ├── rep1/
│   │   │   ├── md.xtc    # Trajectory coordinates (3.6 GB)
│   │   │   ├── md.tpr    # Topology and parameters (16 MB)
│   │   │   ├── md.edr    # Energy data (1.4 MB)
│   │   │   ├── md.log    # Simulation log (1.3 MB)
│   │   │   └── md.gro    # Final structure (35 MB)
│   │   ├── rep2/
│   │   └── rep3/
│   ├── G421R/
│   │   └── [rep1-3]
│   └── R138G/
│       └── [rep1-3]
│
├── benign/
│   └── I300V/
│       └── [rep1-3]
│
└── wildtype/
    └── WT/
        └── [rep1-3]
```

**Total:** 15 replicates × 5 files = 75 files, ~55 GB

## Download Options

### Option 1: Full Dataset (55 GB)
Download all trajectories for comprehensive analysis.

```bash
# Using wget (Zenodo)
wget https://zenodo.org/record/XXXXXX/files/FBXO11_MD_trajectories.tar.gz
tar -xzf FBXO11_MD_trajectories.tar.gz

# Using gdown (Google Drive)
pip install gdown
gdown https://drive.google.com/uc?id=XXXXXX
tar -xzf FBXO11_MD_trajectories.tar.gz
```

### Option 2: Selective Download
Download only specific variants or replicates.

**Example: D910G only (~11 GB)**
```bash
wget https://zenodo.org/record/XXXXXX/files/D910G_trajectories.tar.gz
tar -xzf D910G_trajectories.tar.gz
```

### Option 3: Summary Data Only (Included in GitHub)
If you don't need full trajectories, use the summary CSV files in this repository:
- `trajectory_summary.csv` - RMSD, RMSF, SASA metrics
- `md_analysis_results.csv` - Comprehensive analysis

## File Descriptions

### md.xtc (3.6 GB per replicate)
**GROMACS compressed trajectory file**
- Contains: Atom coordinates for all frames
- Format: Binary (XTC format)
- Frames: 2001 frames (20 ns at 10 ps intervals)
- Size: ~3.6 GB (compressed)

**Usage:**
```bash
# View trajectory in VMD
vmd md.tpr md.xtc

# Calculate RMSD
echo "4 4" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
```

### md.tpr (16 MB per replicate)
**GROMACS portable binary run input file**
- Contains: Topology, parameters, initial structure
- Format: Binary (TPR format)
- Required for: All GROMACS analysis commands
- Size: ~16 MB

**Usage:**
```bash
# Extract structure
gmx dump -s md.tpr > topology.txt
```

### md.edr (1.4 MB per replicate)
**GROMACS energy file**
- Contains: Energy terms (potential, kinetic, temperature, pressure)
- Format: Binary (EDR format)
- Frequency: 10 ps intervals
- Size: ~1.4 MB

**Usage:**
```bash
# Extract potential energy
echo "Potential" | gmx energy -f md.edr -o potential.xvg
```

### md.log (1.3 MB per replicate)
**GROMACS simulation log**
- Contains: Simulation progress, performance metrics
- Format: Plain text
- Size: ~1.3 MB

**Usage:**
```bash
# Check simulation completion
grep "Performance:" md.log
```

### md.gro (35 MB per replicate)
**GROMACS coordinate file (final structure)**
- Contains: Final atomic coordinates
- Format: Plain text (GRO format)
- Size: ~35 MB

**Usage:**
```bash
# View final structure in PyMOL
pymol md.gro
```

## Verification

### File Integrity Check
MD5 checksums provided for verification:

```bash
md5sum -c md5sums.txt
```

**Expected output:**
```
D910G/rep1/md.xtc: OK
D910G/rep1/md.tpr: OK
...
```

### Completeness Check
Verify all expected files present:

```bash
# Count trajectory files
find trajectories/ -name "*.xtc" | wc -l
# Expected: 15

# Check total size
du -sh trajectories/
# Expected: ~55 GB
```

## System Requirements

**Disk space:** 55 GB minimum (60 GB recommended with working space)
**RAM:** 8 GB minimum for analysis (16 GB recommended)
**Software:** GROMACS 2024.x or later, PyMOL/VMD/ChimeraX (optional)

## Analysis Workflow

Once downloaded, place trajectories in the repository structure:

```bash
cd FBXO11_analysis
mkdir -p data/md_trajectories/trajectories
tar -xzf FBXO11_MD_trajectories.tar.gz -C data/md_trajectories/trajectories/

# Verify structure
ls data/md_trajectories/trajectories/pathogenic/D910G/rep1/
# Expected: md.xtc md.tpr md.edr md.log md.gro
```

Then run analysis scripts:

```bash
cd scripts/md_simulations
python 03_analyze_trajectories.py --input ../../data/md_trajectories/trajectories/
```

## Citation

If you use these trajectories, please cite:

**Manuscript:**
Shim Y (2025). Structural Destabilization as a Pathogenic Mechanism in FBXO11-Related Neurodevelopmental Disorder: A Computational Framework for Variant Interpretation. *Frontiers in Genetics*, TBD.

**Dataset:**
Shim Y (2025). FBXO11 Molecular Dynamics Trajectories. Zenodo. [DOI to be assigned]

## Questions

For trajectory access issues or questions:
- Open a GitHub issue: [https://github.com/ykshim2013/FBXO11_analysis/issues](https://github.com/ykshim2013/FBXO11_analysis/issues)
- Contact: ykshim2013@gmail.com

---

**Last updated**: 2025-11-14
