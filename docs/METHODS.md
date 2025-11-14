# Computational Methods

Complete protocols for FBXO11 variant stability analysis.

## Table of Contents

1. [Variant Dataset Assembly](#variant-dataset-assembly)
2. [AlphaFold3 Structure Prediction](#alphafold3-structure-prediction)
3. [FoldX Stability Predictions](#foldx-stability-predictions)
4. [Rosetta Stability Predictions](#rosetta-stability-predictions)
5. [Molecular Dynamics Simulations](#molecular-dynamics-simulations)
6. [Statistical Analysis](#statistical-analysis)

---

## Variant Dataset Assembly

### Pathogenic Variants (n=24)

Pathogenic variants were curated from:
- Published literature (PubMed, Genetics databases)
- ClinVar (Pathogenic or Likely Pathogenic classification)
- DECIPHER database
- HGMD (Human Gene Mutation Database)

**Inclusion criteria:**
- Missense variants (single amino acid substitutions)
- Clinical diagnosis of FBXO11-related neurodevelopmental disorder
- Segregation data or de novo occurrence when available

**Note:** FBXO11-related disorder is a rare Mendelian condition with approximately 50-60 affected individuals reported to date. Our dataset represents near-comprehensive coverage of all clinically validated missense variants.

### Benign Variants (n=21)

Benign controls selected to match pathogenic variants by:
- Domain distribution
- Amino acid chemical properties
- Evolutionary conservation scores
- Surface accessibility

**Sources:**
- gnomAD v4 (high-frequency variants, MAF > 0.1%)
- ClinVar (Benign or Likely Benign when available)
- Conservative substitutions at equivalent positions

---

## AlphaFold3 Structure Prediction

### Server Configuration

**Platform:** AlphaFold3 Server (https://alphafoldserver.com)
**Version:** AlphaFold3 (May 2024 release)
**Input:** FBXO11 amino acid sequence (UniProt: Q86XK2)

### Prediction Settings

**Monomer structures:**
- Sequence: FBXO11 full-length (943 residues)
- Models generated: 5 (ranked by pLDDT)
- Selected seeds: rank_0, rank_3, rank_4

**Complex structures (FBXO11-SKP1):**
- Chain A: FBXO11 (943 residues)
- Chain B: SKP1 (163 residues, UniProt: P63208)
- Models generated: 5 (ranked by interface quality)
- Selected seeds: rank_001, rank_003, rank_005

### Conformational Seed Selection

Three conformational seeds selected per biological context based on:
1. **High structural quality** (pLDDT > 70 for F-box domain)
2. **Structural diversity** (pairwise RMSD > 19 Å)
3. **Interface confidence** (pLDDT > 80 for FBXO11-SKP1 interface)

**Rationale:** Multi-conformational analysis captures structural heterogeneity and reduces overfitting to single structural snapshot.

---

## FoldX Stability Predictions

### Software

**Version:** FoldX 5.1
**Force field:** Empirical energy function (Guerois et al., 2002)
**Platform:** Linux x86_64

### Protocol

**1. Structure Preparation**
```bash
foldx --command=RepairPDB --pdb=structure.pdb
```
- Optimizes rotamers
- Fixes steric clashes
- Standardizes residue geometry

**2. ΔΔG Calculation**
```bash
foldx --command=BuildModel --pdb=structure_Repair.pdb --mutant-file=individual_list.txt
```

**Mutation file format:**
```
K135R;
A486P;
...
```

**3. Technical Replicates**
- 10 replicates per system using different random seeds
- Each replicate performs independent rotamer optimization

**4. ΔΔG Extraction**
```bash
ΔΔG = ΔG_mutant - ΔG_wildtype
```

### Quality Control

- Discard replicates with steric clashes (FoldX warning flags)
- Verify energy convergence (SD < 2.0 kcal/mol across replicates)
- Manual inspection of high-magnitude predictions (|ΔΔG| > 10 kcal/mol)

### Interpretation Thresholds

- **ΔΔG > 1.0 kcal/mol:** Destabilizing (likely pathogenic)
- **-1.0 ≤ ΔΔG ≤ 1.0:** Neutral (VUS)
- **ΔΔG < -1.0:** Stabilizing (likely benign)

---

## Rosetta Stability Predictions

### Software

**Version:** Rosetta 2023.49
**Protocol:** cartesian_ddg
**Force field:** REF2015 (Alford et al., 2017)
**Platform:** Linux x86_64

### Protocol

**1. Structure Relaxation**
```bash
rosetta_scripts.default.linuxgccrelease \
  -s structure.pdb \
  -parser:protocol relax.xml \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:ramp_constraints false
```

**2. Cartesian ΔΔG Calculation**
```bash
cartesian_ddg.default.linuxgccrelease \
  -s structure_relaxed.pdb \
  -ddg::mut_file mutations.txt \
  -ddg::iterations 10 \
  -ddg::cartesian \
  -ddg::score_cutoff 1.0
```

**Mutation file format:**
```
total 1
1
K135 R
```

**3. Technical Replicates**
- 10 replicates per system
- Each replicate uses independent starting conformation from relaxation

**4. ΔΔG Calculation**
```
ΔΔG = E_mutant - E_wildtype
```

### Quality Control

- Verify Rosetta energy convergence (score SD < 5.0 REU)
- Exclude replicates with unrealistic conformations (clash detection)
- Manual inspection of extreme predictions (|ΔΔG| > 20 kcal/mol)

### Interpretation Thresholds

- **ΔΔG > 2.0 kcal/mol:** Strongly destabilizing (likely pathogenic)
- **1.0 < ΔΔG ≤ 2.0:** Moderately destabilizing (VUS/pathogenic)
- **-1.0 ≤ ΔΔG ≤ 1.0:** Neutral (VUS)
- **ΔΔG < -1.0:** Stabilizing (likely benign)

---

## Molecular Dynamics Simulations

### Software

**Version:** GROMACS 2024.4
**Hardware:** NVIDIA GPU (CUDA acceleration)
**Force field:** AMBER99SB-ILDN
**Water model:** TIP3P

### Variant Selection

5 representative variants selected:
- **Pathogenic (n=3):** D910G, G421R, R138G
  - Representing different predicted destabilization severities
- **Benign (n=1):** I300V
  - Conservative substitution control
- **Wildtype (n=1):** WT
  - Reference baseline

### System Preparation

**1. Structure Processing**
```bash
gmx pdb2gmx -f structure.pdb -o processed.gro -water tip3p -ff amber99sb-ildn
```

**2. Solvation**
```bash
gmx editconf -f processed.gro -o box.gro -c -d 1.0 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
```

**3. Ion Addition (Neutralization)**
```bash
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral
```

### Energy Minimization

```bash
gmx grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

**Parameters:**
- Algorithm: Steepest descent
- Max steps: 50,000
- Energy tolerance: 1000.0 kJ/mol/nm

### Equilibration

**NVT Equilibration (100 ps):**
```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
```
- Temperature: 310 K (V-rescale thermostat)
- Coupling constant: 0.1 ps

**NPT Equilibration (100 ps):**
```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
```
- Pressure: 1 bar (Parrinello-Rahman barostat)
- Coupling constant: 2.0 ps

### Production MD

```bash
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -deffnm md -nb gpu
```

**Parameters:**
- Duration: 20 ns per replicate
- Timestep: 2 fs
- Output frequency: 10 ps (2001 frames)
- Temperature: 310 K (V-rescale)
- Pressure: 1 bar (Parrinello-Rahman)
- Constraints: LINCS (bonds involving hydrogen)
- Electrostatics: PME (1.0 nm cutoff)
- Van der Waals: 1.0 nm cutoff

**Replicates:** 3 independent replicates per variant (different initial velocities)

### Trajectory Analysis

**RMSD (Backbone stability):**
```bash
echo "4 4" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns
```

**RMSF (Per-residue flexibility):**
```bash
echo "1" | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res
```

**SASA (Surface accessibility):**
```bash
echo "1" | gmx sasa -s md.tpr -f md.xtc -o sasa.xvg
```

**Radius of gyration (Compactness):**
```bash
echo "1" | gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
```

---

## Statistical Analysis

### Parametric Tests

**Welch's t-test:**
- Comparing mean ΔΔG between pathogenic and benign
- Assumes unequal variances
- Two-tailed test, α = 0.05

**Pearson correlation:**
- FoldX-Rosetta correlation analysis
- Assesses agreement between methods

### Non-Parametric Tests

**Mann-Whitney U test:**
- Comparing median ΔΔG distributions
- Robust to outliers
- One-tailed test (pathogenic > benign), α = 0.05

### Effect Size

**Cohen's d:**
```
d = (mean_pathogenic - mean_benign) / pooled_SD
```

Interpretation:
- |d| < 0.2: Small effect
- 0.2 ≤ |d| < 0.8: Medium effect
- |d| ≥ 0.8: Large effect

### Classification Performance

**ROC/PR Curves:**
- Threshold optimization for ΔΔG cutoffs
- AUC calculation (area under curve)
- Precision-recall analysis

**Combined Model:**
- Logistic regression with FoldX + Rosetta as features
- 10-fold cross-validation
- Performance metrics: AUC, sensitivity, specificity

---

## Software Citations

**AlphaFold3:**
Abramson J, et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*, 630:493-500.

**FoldX:**
Delgado J, et al. (2019). FoldX 5.0: working with RNA, small molecules and a new graphical interface. *Bioinformatics*, 35(20):4168-4169.

Guerois R, Nielsen JE, Serrano L (2002). Predicting changes in the stability of proteins and protein complexes: a study of more than 1000 mutations. *J Mol Biol*, 320(2):369-387.

**Rosetta:**
Alford RF, et al. (2017). The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design. *J Chem Theory Comput*, 13(6):3031-3048.

Park H, et al. (2016). Simultaneous Optimization of Biomolecular Energy Functions on Features from Small Molecules and Macromolecules. *J Chem Theory Comput*, 12(12):6201-6212.

**GROMACS:**
Abraham MJ, et al. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX*, 1:19-25.

**AMBER Force Field:**
Lindorff-Larsen K, et al. (2010). Improved side-chain torsion potentials for the Amber ff99SB protein force field. *Proteins*, 78(8):1950-8.

---

**Last updated**: 2025-11-14
