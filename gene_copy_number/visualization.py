import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import re

class CopyNumberVisualizer:
    def __init__(self):
        pass

    def plot_genome_coverage(self, coverage_file, gene_positions_file, output_file,
                             expected_ploidy=2.0, highlight_genes=None, alpha=0.05):

        highlight_genes = highlight_genes or []
        if any(isinstance(g, list) for g in highlight_genes):
            highlight_genes = [item for sublist in highlight_genes for item in (sublist if isinstance(sublist, list) else [sublist])]

        def load_coverage_data(path):
            df = pd.read_csv(path)
            df = df.rename(columns={  # Normalize mean column
            col: "mean_coverage" for col in df.columns if "mean_cvg" in col.lower()
            })

            # Split 'Target' into chr, start, end
            target_split = df["Target"].str.extract(r'(?P<chr>[^:]+):(?P<start>\d+)-(?P<end>\d+)')
            df = pd.concat([df, target_split], axis=1)
            df[["start", "end"]] = df[["start", "end"]].astype(int)

            return df
        
        def load_gene_positions(path):
            df = pd.read_csv(path, sep='\t', names=['chr', 'start', 'end', 'gene_name', 'gene_id'])
            df['gene_name'] = df['gene_name'].fillna(df['gene_id'])  # fallback to ID if Name is missing
            return df
        
        def calculate_copy_number(df, ploidy):
            # Step 1: Normalize mean coverage to copy number
            median_cov = np.median(df["mean_coverage"])
            df["estimated_cn"] = (df["mean_coverage"] / median_cov) * ploidy

            # Step 2: Remove outliers using IQR approach
            Q1 = df["estimated_cn"].quantile(0.25)
            Q3 = df["estimated_cn"].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            higher_bound = Q3 + 1.5 * IQR

            # Filter out the outliers
            filtered_df = df[(df["estimated_cn"] >= lower_bound) & (df["estimated_cn"] <= higher_bound)]

            # Step 3: Calculate standard z-scores on the filtered_data
            mean_cn = np.mean(filtered_df["estimated_cn"])
            std_cn = np.std(filtered_df["estimated_cn"])
            filtered_df["cn_zscore"] = (filtered_df["estimated_cn"] - mean_cn) / std_cn if std_cn > 0 else 0

            # Step 4: Calculate two-sided p-values from z-scores and an calculate FDR
            filtered_df["p_value"] = 2 * (1 - norm.cdf(np.abs(filtered_df["cn_zscore"])))
            filtered_df["p_adj"] = multipletests(filtered_df["p_value"], method="fdr_bh")[1]
            
            return filtered_df
        
        def merge_data(coverage_df, position_df):
            return pd.merge(coverage_df, position_df, on=["chr", "start", "end"], how="inner")
        
        def natural_sort_key(s):
            return [int(text) if text.isdigit() else text.lower() for text in re.split(r'([0-9]+)', s)]
        
        def roman_sort_key(chrom):
            roman_to_int = {
            'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8,
            'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16
            }
            return roman_to_int.get(chrom.replace('chr', '').upper(), 999)
        
        cov_df = load_coverage_data(coverage_file)
        pos_df = load_gene_positions(gene_positions_file)
        cov_df = calculate_copy_number(cov_df, expected_ploidy)
        merged_df = merge_data(cov_df, pos_df)

        merged_df["chr"] = merged_df["chr"].astype(str)
        example_chr = merged_df["chr"].dropna().iloc[0]
        if re.fullmatch(r'(chr)?[IVX]+', example_chr.upper()):
            merged_df['chr_sort_key'] = merged_df['chr'].map(roman_sort_key)
        else:
            merged_df['chr_sort_key'] = merged_df['chr'].map(natural_sort_key)
        merged_df = merged_df.sort_values(by=['chr_sort_key', 'start']).drop(columns='chr_sort_key')
        merged_df['x'] = range(len(merged_df))
        
        mean_cn = np.mean(merged_df["estimated_cn"])
        std_cn = np.std(merged_df["estimated_cn"])
        # Convert alpha to equivalent z-score
        z_threshold_equivalent = norm.ppf(1 - alpha / 2)
        threshold_high = mean_cn + z_threshold_equivalent * std_cn
        threshold_low = mean_cn - z_threshold_equivalent * std_cn

        # Identify significant genes
        sig_df_raw = merged_df[merged_df["p_value"] < alpha]
        sig_df_fdr = merged_df[merged_df["p_adj"] < alpha]
        sig_df = pd.concat([sig_df_raw, sig_df_fdr]).drop_duplicates()

        print(f"Genes with raw p < 0.05: {(merged_df['p_value'] < 0.05).sum()}")
        print(f"Genes with FDR-adjusted p < {alpha}: {(merged_df['p_adj'] < alpha).sum()}")

        plt.figure(figsize=(16, 6))
        plt.suptitle("Genome-wide Copy Number Variation", fontsize=14)
        sns.scatterplot(data=merged_df, x="x", y="estimated_cn",
                        hue=merged_df["chr"], palette="tab20", legend=False, s=10)
        
        plt.axhline(threshold_high, color='gray', linestyle='--', linewidth=1.0)
        plt.axhline(threshold_low, color='gray', linestyle='--', linewidth=1.0)
        
        # Add chromosome labels and vertical lines between chromosomes
        chrom_groups = merged_df.groupby("chr")["x"]
        y_chrom_label = 1.001 * plt.ylim()[1]  # dynamic top of y-axis - padding
        for chrom, positions in chrom_groups:
            xmin = positions.min()
            xmax = positions.max()
            midpoint = (xmin + xmax) / 2
            plt.axvline(x=xmin, color="gray", linestyle="dotted", linewidth=1.0)
            plt.text(midpoint, y_chrom_label, chrom, ha="center", va="bottom", fontsize=9, fontweight="bold")
        
        # Highlight selected genes
        if highlight_genes:
            highlighted = merged_df[merged_df["gene_name"].isin(highlight_genes)]
            plt.scatter(highlighted["x"], highlighted["estimated_cn"], color="red",
                        s=15, label="Highlighted genes")
            for _, row in highlighted.iterrows():
                plt.text(row["x"], row["estimated_cn"] + 0.2, row["gene_name"], fontsize=8, ha="center")
        
        # Plot significant genes
        for _, row in sig_df.iterrows():
            if row["p_adj"] < alpha and row["gene_name"] in highlight_genes:
                color = "purple"
                label = "Significant & highlighted"
            elif row["p_adj"] < alpha:
                color = "blue"
                label = f"FDR < α = {alpha}"
            elif row["p_value"] < alpha:
                color = "gray"
                label = f"p < α = {alpha}"
            else:
                continue
            offset = 0.05 if row["estimated_cn"] >= expected_ploidy else -0.05
            va = "bottom" if offset > 0 else "top"
            plt.scatter(row["x"], row["estimated_cn"], color=color, s=15, label=label)
            plt.text(row["x"], row["estimated_cn"] + offset, row["gene_name"],
                     fontsize=7, ha="center", va=va, color=color)
            
        custom_lines = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=6, label=f'FDR < α={alpha}'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=6, label=f'p < α={alpha}'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=6, label='Significant & Highlighted')
        ]
        plt.legend(handles=custom_lines, loc='center right', frameon=True)
            
        plt.axhline(expected_ploidy, color="black", linestyle="--", lw=1)
        plt.ylabel("Estimated Copy Number")
        plt.xticks([], [])
        plt.xlabel("")
        ax = plt.gca()
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()

        sig_output = os.path.splitext(output_file)[0] + "_significant_genes.tsv"
        sig_df[["gene_name", "estimated_cn", "cn_zscore", "p_value", "p_adj", "chr", "start", "end"]].to_csv(sig_output, sep="\t", index=False)

        print(f"Plot saved to: {output_file}")
        print(f"Significant genes saved to: {sig_output}")