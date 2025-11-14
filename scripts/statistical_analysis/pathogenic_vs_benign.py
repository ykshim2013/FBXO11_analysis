#!/usr/bin/env python3
"""
Statistical comparison of pathogenic vs benign FBXO11 variants.

Performs Mann-Whitney U test and Welch's t-test to compare ΔΔG distributions
between pathogenic and benign variants for both FoldX and Rosetta predictions.

Author: Youngkyu Shim
Date: 2025-11-14
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(foldx_path, rosetta_path):
    """Load FoldX and Rosetta aggregated results."""
    foldx_df = pd.read_csv(foldx_path)
    rosetta_df = pd.read_csv(rosetta_path)
    return foldx_df, rosetta_df

def mann_whitney_test(pathogenic, benign, method_name):
    """Perform Mann-Whitney U test between pathogenic and benign."""
    statistic, p_value = stats.mannwhitneyu(pathogenic, benign, alternative='greater')

    # Calculate effect size (Cohen's d)
    n1, n2 = len(pathogenic), len(benign)
    var1, var2 = np.var(pathogenic, ddof=1), np.var(benign, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    cohens_d = (np.mean(pathogenic) - np.mean(benign)) / pooled_std

    print(f"\n{method_name} - Mann-Whitney U Test:")
    print(f"  Pathogenic median: {np.median(pathogenic):.2f} kcal/mol")
    print(f"  Benign median: {np.median(benign):.2f} kcal/mol")
    print(f"  U-statistic: {statistic:.2f}")
    print(f"  p-value: {p_value:.2e}")
    print(f"  Cohen's d: {cohens_d:.2f}")

    return {'statistic': statistic, 'p_value': p_value, 'cohens_d': cohens_d}

def welch_ttest(pathogenic, benign, method_name):
    """Perform Welch's t-test between pathogenic and benign."""
    statistic, p_value = stats.ttest_ind(pathogenic, benign, equal_var=False)

    # Calculate effect size (Cohen's d)
    n1, n2 = len(pathogenic), len(benign)
    var1, var2 = np.var(pathogenic, ddof=1), np.var(benign, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    cohens_d = (np.mean(pathogenic) - np.mean(benign)) / pooled_std

    print(f"\n{method_name} - Welch's t-test:")
    print(f"  Pathogenic mean: {np.mean(pathogenic):.2f} ± {np.std(pathogenic):.2f} kcal/mol")
    print(f"  Benign mean: {np.mean(benign):.2f} ± {np.std(benign):.2f} kcal/mol")
    print(f"  t-statistic: {statistic:.2f}")
    print(f"  p-value: {p_value:.2e}")
    print(f"  Cohen's d: {cohens_d:.2f}")

    return {'statistic': statistic, 'p_value': p_value, 'cohens_d': cohens_d}

def plot_distributions(foldx_df, rosetta_df, output_path):
    """Create violin plots comparing pathogenic vs benign distributions."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # FoldX plot
    sns.violinplot(data=foldx_df, x='classification', y='ddg', ax=axes[0],
                   order=['pathogenic', 'benign'], palette=['#e74c3c', '#3498db'])
    axes[0].set_title('FoldX ΔΔG Distribution', fontsize=14, fontweight='bold')
    axes[0].set_xlabel('Classification', fontsize=12)
    axes[0].set_ylabel('ΔΔG (kcal/mol)', fontsize=12)
    axes[0].axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Threshold')
    axes[0].legend()

    # Rosetta plot
    sns.violinplot(data=rosetta_df, x='classification', y='ddg', ax=axes[1],
                   order=['pathogenic', 'benign'], palette=['#e74c3c', '#3498db'])
    axes[1].set_title('Rosetta ΔΔG Distribution', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('Classification', fontsize=12)
    axes[1].set_ylabel('ΔΔG (kcal/mol)', fontsize=12)
    axes[1].axhline(y=2.0, color='gray', linestyle='--', alpha=0.5, label='Threshold')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nDistribution plot saved to: {output_path}")

def main():
    """Main analysis workflow."""
    print("="*60)
    print("FBXO11 Pathogenic vs Benign Statistical Analysis")
    print("="*60)

    # Load data
    foldx_df, rosetta_df = load_data(
        '../../data/foldx/aggregated_results.csv',
        '../../data/rosetta/aggregated_results.csv'
    )

    # Filter pathogenic and benign (exclude wildtype)
    foldx_path = foldx_df[foldx_df['classification'] == 'pathogenic']['ddg']
    foldx_benign = foldx_df[foldx_df['classification'] == 'benign']['ddg']

    rosetta_path = rosetta_df[rosetta_df['classification'] == 'pathogenic']['ddg']
    rosetta_benign = rosetta_df[rosetta_df['classification'] == 'benign']['ddg']

    # FoldX analysis
    foldx_mw = mann_whitney_test(foldx_path, foldx_benign, "FoldX")
    foldx_welch = welch_ttest(foldx_path, foldx_benign, "FoldX")

    # Rosetta analysis
    rosetta_mw = mann_whitney_test(rosetta_path, rosetta_benign, "Rosetta")
    rosetta_welch = welch_ttest(rosetta_path, rosetta_benign, "Rosetta")

    # Create visualizations
    plot_distributions(foldx_df, rosetta_df, '../../analysis_outputs/figures/pathogenic_vs_benign.png')

    # Summary table
    print("\n" + "="*60)
    print("Summary Table")
    print("="*60)
    print(f"{'Method':<10} {'Test':<15} {'p-value':<12} {'Cohen's d':<10}")
    print("-"*60)
    print(f"{'FoldX':<10} {'Mann-Whitney':<15} {foldx_mw['p_value']:<12.2e} {foldx_mw['cohens_d']:<10.2f}")
    print(f"{'FoldX':<10} {'Welch t-test':<15} {foldx_welch['p_value']:<12.2e} {foldx_welch['cohens_d']:<10.2f}")
    print(f"{'Rosetta':<10} {'Mann-Whitney':<15} {rosetta_mw['p_value']:<12.2e} {rosetta_mw['cohens_d']:<10.2f}")
    print(f"{'Rosetta':<10} {'Welch t-test':<15} {rosetta_welch['p_value']:<12.2e} {rosetta_welch['cohens_d']:<10.2f}")
    print("="*60)

if __name__ == '__main__':
    main()
