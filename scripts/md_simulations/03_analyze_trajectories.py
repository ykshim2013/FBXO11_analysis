#!/usr/bin/env python3
"""
Molecular dynamics trajectory analysis for FBXO11 variants.

Calculates RMSD, RMSF, SASA, and radius of gyration for all MD trajectories
using MDAnalysis and GROMACS tools.

Author: Youngkyu Shim
Date: 2025-11-14
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, rmsf
import argparse

def run_gromacs_command(command, input_text=None):
    """Run GROMACS command with optional piped input."""
    try:
        if input_text:
            result = subprocess.run(
                command,
                input=input_text,
                capture_output=True,
                text=True,
                check=True
            )
        else:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True
            )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(command)}")
        print(f"Error message: {e.stderr}")
        sys.exit(1)

def calculate_rmsd(tpr_file, xtc_file, output_file):
    """Calculate backbone RMSD using GROMACS."""
    print(f"  Calculating RMSD...")
    cmd = ['gmx', 'rms', '-s', tpr_file, '-f', xtc_file, '-o', output_file, '-tu', 'ns']
    run_gromacs_command(cmd, input_text="4 4\n")  # Backbone selection

    # Parse RMSD data
    rmsd_data = np.loadtxt(output_file, comments=['@', '#'])
    return {
        'rmsd_mean': np.mean(rmsd_data[:, 1]),
        'rmsd_sd': np.std(rmsd_data[:, 1]),
        'rmsd_max': np.max(rmsd_data[:, 1])
    }

def calculate_rmsf(tpr_file, xtc_file, output_file):
    """Calculate per-residue RMSF using GROMACS."""
    print(f"  Calculating RMSF...")
    cmd = ['gmx', 'rmsf', '-s', tpr_file, '-f', xtc_file, '-o', output_file, '-res']
    run_gromacs_command(cmd, input_text="1\n")  # Protein selection

    # Parse RMSF data
    rmsf_data = np.loadtxt(output_file, comments=['@', '#'])
    return {
        'rmsf_mean': np.mean(rmsf_data[:, 1])
    }

def calculate_sasa(tpr_file, xtc_file, output_file):
    """Calculate solvent-accessible surface area using GROMACS."""
    print(f"  Calculating SASA...")
    cmd = ['gmx', 'sasa', '-s', tpr_file, '-f', xtc_file, '-o', output_file]
    run_gromacs_command(cmd, input_text="1\n")  # Protein selection

    # Parse SASA data
    sasa_data = np.loadtxt(output_file, comments=['@', '#'])
    return {
        'sasa_mean': np.mean(sasa_data[:, 1])
    }

def calculate_rg(tpr_file, xtc_file, output_file):
    """Calculate radius of gyration using GROMACS."""
    print(f"  Calculating radius of gyration...")
    cmd = ['gmx', 'gyrate', '-s', tpr_file, '-f', xtc_file, '-o', output_file]
    run_gromacs_command(cmd, input_text="1\n")  # Protein selection

    # Parse Rg data
    rg_data = np.loadtxt(output_file, comments=['@', '#'])
    return {
        'rg_mean': np.mean(rg_data[:, 1])
    }

def analyze_trajectory(trajectory_dir, variant, classification, replicate):
    """Analyze a single MD trajectory."""
    print(f"\nAnalyzing {variant} replicate {replicate}...")

    tpr_file = os.path.join(trajectory_dir, 'md.tpr')
    xtc_file = os.path.join(trajectory_dir, 'md.xtc')

    # Check if files exist
    if not os.path.exists(tpr_file) or not os.path.exists(xtc_file):
        print(f"  Warning: Missing trajectory files in {trajectory_dir}")
        return None

    # Create output directory
    output_dir = os.path.join(trajectory_dir, 'analysis')
    os.makedirs(output_dir, exist_ok=True)

    # Calculate metrics
    results = {
        'variant': variant,
        'classification': classification,
        'replicate': replicate
    }

    try:
        rmsd_results = calculate_rmsd(tpr_file, xtc_file, os.path.join(output_dir, 'rmsd.xvg'))
        rmsf_results = calculate_rmsf(tpr_file, xtc_file, os.path.join(output_dir, 'rmsf.xvg'))
        sasa_results = calculate_sasa(tpr_file, xtc_file, os.path.join(output_dir, 'sasa.xvg'))
        rg_results = calculate_rg(tpr_file, xtc_file, os.path.join(output_dir, 'gyrate.xvg'))

        results.update(rmsd_results)
        results.update(rmsf_results)
        results.update(sasa_results)
        results.update(rg_results)

    except Exception as e:
        print(f"  Error analyzing trajectory: {e}")
        return None

    print(f"  ✓ Analysis complete")
    return results

def main():
    """Main analysis workflow."""
    parser = argparse.ArgumentParser(description='Analyze FBXO11 MD trajectories')
    parser.add_argument('--input', required=True, help='Path to trajectories directory')
    parser.add_argument('--output', default='../../data/md_trajectories/', help='Output directory')
    args = parser.parse_args()

    print("="*60)
    print("FBXO11 MD Trajectory Analysis")
    print("="*60)

    # Define trajectory structure
    variants = {
        'pathogenic': ['D910G', 'G421R', 'R138G'],
        'benign': ['I300V'],
        'wildtype': ['WT']
    }

    all_results = []

    # Analyze all trajectories
    for classification, variant_list in variants.items():
        for variant in variant_list:
            for replicate in [1, 2, 3]:
                trajectory_dir = os.path.join(
                    args.input,
                    classification,
                    variant,
                    f'rep{replicate}'
                )

                results = analyze_trajectory(trajectory_dir, variant, classification, replicate)
                if results:
                    all_results.append(results)

    # Create summary dataframe
    if all_results:
        df = pd.DataFrame(all_results)
        output_file = os.path.join(args.output, 'trajectory_summary.csv')
        df.to_csv(output_file, index=False)
        print(f"\n✓ Summary saved to: {output_file}")

        # Print summary statistics
        print("\n" + "="*60)
        print("Summary Statistics")
        print("="*60)
        for classification in df['classification'].unique():
            subset = df[df['classification'] == classification]
            print(f"\n{classification.upper()}:")
            for variant in subset['variant'].unique():
                var_subset = subset[subset['variant'] == variant]
                print(f"  {variant}:")
                print(f"    RMSD: {var_subset['rmsd_mean'].mean():.2f} ± {var_subset['rmsd_mean'].std():.2f} nm")
                print(f"    RMSF: {var_subset['rmsf_mean'].mean():.4f} ± {var_subset['rmsf_mean'].std():.4f} nm")
                print(f"    SASA: {var_subset['sasa_mean'].mean():.2f} ± {var_subset['sasa_mean'].std():.2f} nm²")
    else:
        print("\n✗ No valid trajectories analyzed")

    print("="*60)

if __name__ == '__main__':
    main()
