import pandas as pd
import numpy as np
import glob
import os
import sys
import argparse
from pathlib import Path


def parse_nanostats_file(filepath):
    """Parse a NanoStats.txt file and return metrics as a dictionary"""
    metrics = {}

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if '\t' in line:
                parts = line.split('\t')
                if len(parts) == 2:
                    key, value = parts
                    try:
                        if '.' in value:
                            metrics[key] = float(value)
                        else:
                            metrics[key] = int(value)
                    except ValueError:
                        metrics[key] = value

    return metrics


def load_multiqc_stats(results_dir):
    """Load MultiQC general stats"""
    multiqc_stats_file = results_dir / \
        "multiqc/minimap2/multiqc_data/multiqc_general_stats.txt"

    if not multiqc_stats_file.exists():
        print(f"MultiQC stats file not found: {multiqc_stats_file}")
        return pd.DataFrame()

    multiqc_df = pd.read_csv(multiqc_stats_file, sep='\t', index_col=0)
    print(
        f"Loaded MultiQC stats: {multiqc_df.shape[0]} samples, {multiqc_df.shape[1]} metrics")
    return multiqc_df


def load_nanoplot_metrics(results_dir):
    """Extract NanoPlot metrics from NanoStats.txt files"""
    nanostats_files = list(results_dir.glob("nanoplot/*NanoStats.txt"))
    print(f"Found {len(nanostats_files)} NanoStats files")

    if not nanostats_files:
        print("Warning: No NanoStats files found")
        return pd.DataFrame()

    nanoplot_metrics = {}
    for file in nanostats_files:
        sample_name = file.name.replace(
            'NanoStats.txt', '').replace('_post_filtering', '')
        metrics = parse_nanostats_file(file)
        nanoplot_metrics[sample_name] = metrics

        print(f"  Parsed metrics for: {sample_name}")

    nanoplot_df = pd.DataFrame(nanoplot_metrics).T
    print(
        f"NanoPlot metrics: {nanoplot_df.shape[0]} samples, {nanoplot_df.shape[1]} metrics")
    return nanoplot_df


def calculate_gene_metrics(results_dir):
    """Calculate gene expression summary metrics from featureCounts"""
    # Try bambu first, then stringtie2
    gene_counts_file = results_dir / "bambu/counts_gene.txt"
    if not gene_counts_file.exists():
        gene_counts_file = results_dir / "stringtie2/featureCounts/counts_gene.txt"

    if not gene_counts_file.exists():
        print(f"Warning: Gene counts file not found: {gene_counts_file}")
        return pd.DataFrame()

    gene_counts_df = pd.read_csv(gene_counts_file, sep='\t', comment='#')
    print(f"Loaded gene counts: {gene_counts_df.shape[0]} entries")

    gene_id_col = 'GENEID' if 'GENEID' in gene_counts_df.columns else 'Geneid'

    sample_cols = [col for col in gene_counts_df.columns if col.endswith('.sorted.bam')]
    if not sample_cols:
        sample_cols = [col for col in gene_counts_df.columns if col.endswith('.sorted')]

    gene_grouped = gene_counts_df.groupby(gene_id_col)[sample_cols].sum()
    print(f"Collapsed to {gene_grouped.shape[0]} unique genes")

    metrics = {}
    for sample_col in sample_cols:
        sample_name = sample_col.replace('.sorted.bam', '').replace('.sorted', '')

        counts = gene_grouped[sample_col]
        metrics[sample_name] = {
            'genes_detected': int((counts > 0).sum()),
            'genes_above_5_counts': int((counts > 5).sum()),
            'genes_above_10_counts': int((counts > 10).sum()),
            'genes_above_50_counts': int((counts > 50).sum()),
            'total_gene_assigned_reads': int(counts.sum()),
            'median_gene_expression': float(counts[counts > 0].median()) if (counts > 0).any() else 0.0,
            'mean_gene_expression': float(counts[counts > 0].mean()) if (counts > 0).any() else 0.0
        }
        print(
            f"  Gene-level metrics for {sample_name}: {metrics[sample_name]['genes_detected']} genes detected")

    return pd.DataFrame(metrics).T


def combine_metrics(multiqc_df, nanoplot_df, gene_metrics_df):
    """Combine all metrics into a single dataframe"""
    if multiqc_df.empty:
        print("Error: No MultiQC data to use as base")
        return pd.DataFrame()

    combined_df = multiqc_df.copy()
    print(f"Base MultiQC samples: {list(combined_df.index)}")

    if not nanoplot_df.empty:
        nanoplot_cols_to_include = [
            'number_of_reads', 'number_of_bases', 'median_read_length',
            'mean_read_length', 'read_length_stdev', 'n50', 'mean_qual', 'median_qual'
        ]

        added_cols = 0
        for col in nanoplot_cols_to_include:
            if col in nanoplot_df.columns:
                for multiqc_sample in combined_df.index:
                    # Find matching NanoPlot sample
                    matching_nanoplot = None
                    for nanoplot_sample in nanoplot_df.index:
                        if multiqc_sample in nanoplot_sample or nanoplot_sample in multiqc_sample:
                            matching_nanoplot = nanoplot_sample
                            break

                    if matching_nanoplot and col in nanoplot_df.columns:
                        combined_df.loc[multiqc_sample,
                                        f'nanoplot_{col}'] = nanoplot_df.loc[matching_nanoplot, col]
                        added_cols += 1
        print(
            f"Added {len([c for c in combined_df.columns if c.startswith('nanoplot_')])} NanoPlot columns")

    if not gene_metrics_df.empty:
        for col in gene_metrics_df.columns:
            for multiqc_sample in combined_df.index:
                matching_gene = None
                for gene_sample in gene_metrics_df.index:
                    if multiqc_sample in gene_sample or gene_sample in multiqc_sample:
                        matching_gene = gene_sample
                        break

                if matching_gene:
                    combined_df.loc[multiqc_sample,
                                    f'{col}'] = gene_metrics_df.loc[matching_gene, col]
        print(
            f"Added {len([c for c in combined_df.columns if c.startswith('gene_')])} gene metric columns")

    print(
        f"Combined metrics: {combined_df.shape[0]} samples, {combined_df.shape[1]} metrics")
    return combined_df




def clean_dataframe_for_export(df, index_name="sample"):
    """Clean dataframe column names for CSV export"""
    df_clean = df.copy()

    clean_columns = {}
    for col in df_clean.columns:
        clean_col = col.replace('(', '_').replace(')', '_').replace(
            ':', '_').replace('>', '_gt_').replace(' ', '_').replace('-', '_')
        clean_col = '_'.join([part for part in clean_col.split('_') if part])
        clean_columns[col] = clean_col

    df_clean = df_clean.rename(columns=clean_columns)

    if df_clean.index.name is None:
        df_clean.index.name = index_name

    return df_clean


def export_results(combined_df, nanoplot_df, gene_metrics_df, results_dir):
    """Export results to CSV files"""
    output_files = []

    if not combined_df.empty:
        combined_clean = clean_dataframe_for_export(combined_df, "sample")
        output_file = results_dir / "combined_metrics_table.csv"
        combined_clean.to_csv(output_file)
        output_files.append(output_file)
        print(f"Saved combined metrics: {output_file}")

    if not nanoplot_df.empty:
        nanoplot_clean = clean_dataframe_for_export(nanoplot_df, "sample")
        if 'Metrics' in nanoplot_clean.columns:
            nanoplot_clean = nanoplot_clean.drop('Metrics', axis=1)
        output_file = results_dir / "nanoplot_metrics.csv"
        nanoplot_clean.to_csv(output_file)
        output_files.append(output_file)
        print(f"Saved NanoPlot metrics: {output_file}")

    if not gene_metrics_df.empty:
        gene_clean = clean_dataframe_for_export(gene_metrics_df, "sample")
        output_file = results_dir / "gene_expression_metrics.csv"
        gene_clean.to_csv(output_file)
        output_files.append(output_file)
        print(f"Saved gene metrics: {output_file}")

    return output_files


def main():
    parser = argparse.ArgumentParser(
        description='Collate nanoseq metrics into CSV files')
    parser.add_argument('results_dir', help='Results directory')
    parser.add_argument('--output-prefix', default='',
                        help='Prefix for output CSV files')

    args = parser.parse_args()

    results_dir = Path(args.results_dir)

    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        sys.exit(1)

    print(f"\nProcessing results from: {results_dir}\n")

    multiqc_df = load_multiqc_stats(results_dir)
    nanoplot_df = load_nanoplot_metrics(results_dir)
    gene_metrics_df = calculate_gene_metrics(results_dir)

    combined_df = combine_metrics(multiqc_df, nanoplot_df, gene_metrics_df)

    if not combined_df.empty:
        output_files = export_results(combined_df, nanoplot_df, gene_metrics_df, results_dir)
        print(f"\nSuccessfully created {len(output_files)} output file(s)")
        for f in output_files:
            print(f"  - {f}")
    else:
        print("\nError: No metrics to export")
        sys.exit(1)

if __name__ == "__main__":
    main()
