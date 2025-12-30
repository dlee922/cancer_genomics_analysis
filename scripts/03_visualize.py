"""
Visualization Script for TCGA LUAD Analysis
===========================================
Creates publication-quality figures showing mutation patterns and survival analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import os
import sys

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.display import print_header, print_step, print_success, print_info

# Setup
PROCESSED_DATA_DIR = 'data/processed'
OUTPUT_DIR = 'outputs/figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set visualization style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

print_header("TCGA LUAD Visualization")


def load_processed_data():
    """
    Load all processed data files

    Returns: tuple of all processed data files (5)
    """
    print_info("Loading processed data sets")

    print_info("Loading tmb data...")
    tmb_df = pd.read_csv(f'{PROCESSED_DATA_DIR}/tmb_data.csv')

    print_info("Loading mutation matrix data...")
    mutation_matrix = pd.read_csv(f'{PROCESSED_DATA_DIR}/mutation_matrix.csv', index_col=0)
    
    print_info("Loading survival data...")
    survival_df = pd.read_csv(f'{PROCESSED_DATA_DIR}/survival_data.csv')
    
    print_info("Loading integrated data...")
    integrated_df = pd.read_csv(f'{PROCESSED_DATA_DIR}/integrated_data.csv')
    
    print_info("Loading top mutated genes data...")
    top_genes_df = pd.read_csv(f'{PROCESSED_DATA_DIR}/top_mutated_genes.csv')
    
    print_info(f"TMB data: {tmb_df.shape}")
    print_info(f"Mutation matrix data: {mutation_matrix.shape}")
    print_info(f"Survival data: {survival_df.shape}")
    print_info(f"Integrated data: {integrated_df.shape}")
    print_info(f"Top Mutated Genes data: {top_genes_df.shape}")
    
    print_success("All data loaded successfully")
    
    return tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df


def plot_gene_frequency(top_genes_df):
    """
    Create a horizontal bar chart of top mutated genes.

    BIOLOGICAL CONTEXT:
    -------------------
    Shows the "mutation frequency" - what percentage of patients have 
    mutations in each gene. High frequency = likely driver gene.
    
    YOUR TASK:
    ----------
    1. Sort top_genes_df by 'percentage' in descending order
    2. Take the top 15 genes (for visual clarity)
    3. Extract gene names and percentages as lists for plotting
    4. Print what you're plotting (e.g., "Plotting 15 genes...")
    
    Args:
        top_genes_df (pd.DataFrame): Top mutated genes with frequencies
    """
    print_info("Creating gene frequency bar chart...")

    # sort genes by percentage (frequency) and descending order 
    sorted_genes = top_genes_df.sort_values('percentage', ascending=False)

    # grabbing only the top 15 genes so the bar chart is not humungous
    # change based on spefications
    plot_data = sorted_genes.head(15)
    
    # prep data for matplotlib
    genes = plot_data['gene_symbol'].tolist()
    percentages = plot_data['percentage'].tolist()

    print_info(f"Plotting top {len(genes)} mutated genes")

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create horizontal bar chart
    # We reverse the order so highest is at top
    y_positions = np.arange(len(genes))
    bars = ax.barh(y_positions, percentages[::-1], color='steelblue', edgecolor='navy', linewidth=1.2)

    # Customize the plot
    ax.set_yticks(y_positions)
    ax.set_yticklabels(genes[::-1], fontsize=11)  # Reverse to match bars
    ax.set_xlabel('Percentage of Samples Mutated (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {len(genes)} Most Frequently Mutated Genes in LUAD', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add percentage labels on the bars
    for i, (bar, pct) in enumerate(zip(bars, percentages[::-1])):
        # function arguments: (x-coord, y-coord, text to display, etc.)
        ax.text(pct + 1, bar.get_y() + bar.get_height()/2, 
                f'{pct:.1f}%', 
                va='center', fontsize=10, fontweight='bold')
    
    # Add grid for easier reading
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Clean up
    plt.tight_layout()
    
    # Save the figure
    output_path = f'{OUTPUT_DIR}/gene_frequency.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print_success(f"Saved: {output_path}")
    
    # Close to free memory
    plt.close()


def plot_tmb_distribution(tmb_df):
    """
    Create a histogram showing TMB distribution with clinical cutoffs.
    
    BIOLOGICAL CONTEXT:
    -------------------
    TMB is used clinically to predict immunotherapy response:
    - Low TMB (< 10 mut/Mb): Typical response
    - High TMB (≥ 10 mut/Mb): Better response to checkpoint inhibitors
    
    YOUR TASK:
    ----------
    1. Calculate summary statistics: mean, median, std deviation
    2. Define the clinical cutoff (10 mutations/Mb)
    3. Count how many samples are "high TMB" (≥ 10)
    4. Calculate percentage of high TMB samples
    5. Print these statistics
    
    Args:
        tmb_df (pd.DataFrame): TMB data with 'tmb' column
    """
    print_info("Creating TMB distribution plot...")
    
    # HINT: Use .mean(), .median(), .std() on tmb_df['tmb']
    mean_tmb = tmb_df['tmb'].mean()
    median_tmb = tmb_df['tmb'].median()
    std_tmb = tmb_df['tmb'].std()
    
    # clinical cutoff
    TMB_CUTOFF = 10
    
    # TODO: Count high TMB samples
    high_tmb_count = (tmb_df['tmb'] >= TMB_CUTOFF).sum()
    total_samples = len(tmb_df)
    high_tmb_percentage = (high_tmb_count / total_samples) * 100

    print_info(f"Mean TMB: {mean_tmb:.2f} mut/Mb")
    print_info(f"Median TMB: {median_tmb:.2f} mut/Mb")
    print_info(f"High TMB samples (≥{TMB_CUTOFF}): {high_tmb_count} ({high_tmb_percentage:.1f}%)")    
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Create histogram
    # bins=30 means split the data into 30 bins
    # edgecolor makes bars stand out
    # alpha is transparency (0-1)
    n, bins, patches = ax.hist(tmb_df['tmb'], bins=30, 
                                color='skyblue', edgecolor='navy', 
                                alpha=0.7, linewidth=1.2)
    
    # Add vertical line for clinical cutoff
    ax.axvline(TMB_CUTOFF, color='red', linestyle='--', 
               linewidth=2.5, label=f'Clinical Cutoff ({TMB_CUTOFF} mut/Mb)')
    
    # Add vertical lines for mean and median
    ax.axvline(mean_tmb, color='green', linestyle='-', 
               linewidth=2, label=f'Mean ({mean_tmb:.1f} mut/Mb)')
    ax.axvline(median_tmb, color='orange', linestyle='-', 
               linewidth=2, label=f'Median ({median_tmb:.1f} mut/Mb)')
    
    # Add shaded regions for high/low TMB
    ax.axvspan(0, TMB_CUTOFF, alpha=0.1, color='blue', label='Low TMB')
    ax.axvspan(TMB_CUTOFF, tmb_df['tmb'].max(), alpha=0.1, color='red', label='High TMB')
    
    # Labels and title
    ax.set_xlabel('Tumor Mutation Burden (mutations/Mb)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of Samples', fontsize=13, fontweight='bold')
    ax.set_title('Distribution of Tumor Mutation Burden in LUAD\n' + 
                 f'(n={total_samples} samples, {high_tmb_percentage:.1f}% high TMB)',
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add legend
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Add text box with statistics
    stats_text = f'Mean: {mean_tmb:.2f}\nMedian: {median_tmb:.2f}\nStd Dev: {std_tmb:.2f}'
    ax.text(0.02, 0.98, stats_text, 
            transform=ax.transAxes,  # Use axis coordinates (0-1)
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=10, family='monospace')
    
    # Clean up
    plt.tight_layout()
    
    # Save the figure
    output_path = f'{OUTPUT_DIR}/tmb_distribution.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print_success(f"Saved: {output_path}")
    
    plt.close()


def plot_oncoplot(mutation_matrix, top_genes_df):
    """
    Create an oncoplot (heatmap) showing mutation patterns.
    
    BIOLOGICAL CONTEXT:
    -------------------
    The oncoplot is THE signature visualization in cancer genomics.
    - Rows = samples (patients)
    - Columns = genes
    - Blue = mutated, White = wild-type (normal)
    
    Shows patterns like:
    - TP53 + KRAS often mutated together (co-occurrence)
    - EGFR + KRAS rarely both mutated (mutual exclusivity)
    
    YOUR TASK:
    ----------
    1. Select top 15 genes (for visual clarity, 20 is too crowded)
    2. Filter mutation_matrix to only include these genes
    3. Sort samples by total mutation count (most mutated on top)
       HINT: mutation_matrix.sum(axis=1) counts mutations per sample
    4. Make sure genes are ordered by frequency (most mutated on left)
    5. Print matrix dimensions
    
    Args:
        mutation_matrix (pd.DataFrame): Binary mutation matrix (517×20)
        top_genes_df (pd.DataFrame): Gene frequencies
    """
    print_info("Creating oncoplot (mutation heatmap)...")
    
    # select top 15 genes
    top_15_genes = top_genes_df.head(15)['gene_symbol'].tolist()
    
    # filter mutation matrix to these genes
    plot_matrix = mutation_matrix[top_15_genes]
    
    # sort samples by mutation count (most mutated samples on top)
    sample_mutation_counts = plot_matrix.sum(axis=1)  # axis=1 means sum across columns
    sorted_samples = sample_mutation_counts.sort_values(ascending=False).index
    plot_matrix = plot_matrix.loc[sorted_samples]  # Reorder rows
    
    # print info about what we're plotting
    print_info(f"Plotting {plot_matrix.shape[0]} samples × {plot_matrix.shape[1]} genes")
    print_info(f"Genes: {', '.join(top_15_genes[:5])}...")  # Show first 5
    
    # Create figure with subplots
    # We'll have main heatmap + gene frequency bars on top
    fig = plt.figure(figsize=(14, 10))
    
    # Create grid for layout: [gene bars, heatmap]
    # height_ratios controls relative heights
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 10], hspace=0.05)
    
    # Top subplot: gene frequency bars
    ax_genes = fig.add_subplot(gs[0])
    gene_freqs = top_genes_df.head(15)['percentage'].values
    ax_genes.bar(range(len(top_15_genes)), gene_freqs, color='coral', edgecolor='darkred')
    ax_genes.set_xlim(-0.5, len(top_15_genes) - 0.5)
    ax_genes.set_ylabel('Frequency\n(%)', fontsize=10, fontweight='bold')
    ax_genes.set_xticks([])  # No x-axis labels here
    ax_genes.spines['bottom'].set_visible(False)
    ax_genes.grid(axis='y', alpha=0.3)
    
    # Bottom subplot: main heatmap
    ax_heatmap = fig.add_subplot(gs[1])
    
    # Create the heatmap using seaborn
    # cmap: color map (white=0, blue=1)
    # cbar_kws: customize the color bar
    # Create the heatmap using seaborn
    sns.heatmap(plot_matrix, 
                cmap=['#f5f5f5', 'steelblue'],  # Very light gray + blue
                cbar_kws={'label': 'Mutation Status', 
                        'ticks': [0.25, 0.75],
                        'shrink': 0.5},
                linewidths=0.05,  # Very thin lines
                linecolor='#e0e0e0',  # Light gray
                xticklabels=top_15_genes,
                yticklabels=False,
                ax=ax_heatmap)
        
    # Customize color bar labels
    colorbar = ax_heatmap.collections[0].colorbar
    colorbar.ax.set_yticklabels(['Wild-type', 'Mutated'])
    
    # Labels
    ax_heatmap.set_xlabel('Genes', fontsize=12, fontweight='bold')
    ax_heatmap.set_ylabel(f'Samples (n={plot_matrix.shape[0]})', 
                          fontsize=12, fontweight='bold')
    
    # Rotate gene labels
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), 
                               rotation=45, ha='right', fontsize=11)
    
    # Overall title
    fig.suptitle('Oncoplot: Mutation Landscape of Top 15 Genes in LUAD', 
                 fontsize=15, fontweight='bold', y=0.995)
    
    # Add text annotation about sorting
    fig.text(0.99, 0.02, 
             'Samples sorted by mutation count (highest on top)\nBlue = Mutated, White = Wild-type',
             ha='right', fontsize=9, style='italic', color='gray')
    
    # Calculate and display co-mutation statistics
    # Count samples with both TP53 and KRAS mutations (if both in top 15)
    if 'TP53' in top_15_genes and 'KRAS' in top_15_genes:
        both_mutated = ((plot_matrix['TP53'] == 1) & (plot_matrix['KRAS'] == 1)).sum()
        print_info(f"Co-mutation: {both_mutated} samples have both TP53 and KRAS mutations")
    
    # Clean up
    plt.tight_layout()
    
    # Save the figure
    output_path = f'{OUTPUT_DIR}/oncoplot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print_success(f"Saved: {output_path}")
    
    plt.close() 


def plot_survival_curves(survival_df, tmb_df):
    """
    Create Kaplan-Meier survival curves comparing high vs low TMB.
    
    BIOLOGICAL CONTEXT:
    -------------------
    Kaplan-Meier curves show survival probability over time.
    We're asking: "Do patients with high TMB live longer or shorter?"
    
    High TMB might mean:
    - Better immunotherapy response → longer survival
    - OR more aggressive tumor → shorter survival
    
    YOUR TASK:
    ----------
    1. Extract patientId from sampleId in tmb_df (first 12 characters)
    2. Merge survival_df with tmb_df on patientId
    3. Define TMB groups: High (≥10) vs Low (<10)
    4. Count patients in each group
    5. Print group sizes
    
    Args:
        survival_df (pd.DataFrame): Survival data (patientId, time, event)
        tmb_df (pd.DataFrame): TMB data (sampleId, mutation_count, tmb)
    """
    print_info("Creating Kaplan-Meier survival curves...")
    
    # TODO: Extract patientId from sampleId
    # TCGA format: "TCGA-05-4244-01" → patient is "TCGA-05-4244" (first 12 chars)
    # HINT: Use .str[:12] to get first 12 characters
    # tmb_with_patient = tmb_df.copy()
    # tmb_with_patient['patientId'] = tmb_with_patient['sampleId'].str[:12]
    
    # TODO: Merge with survival data
    # HINT: Use pd.merge() or .merge()
    # Keep only patients with both TMB and survival data (inner join)
    # survival_tmb = survival_df.merge(tmb_with_patient[['patientId', 'tmb']], on='patientId', how='inner')
    
    # TODO: Define TMB groups
    # TMB_CUTOFF = 10
    # survival_tmb['tmb_group'] = survival_tmb['tmb'].apply(
    #     lambda x: 'High TMB (≥10)' if x >= TMB_CUTOFF else 'Low TMB (<10)'
    # )
    
    # TODO: Count patients in each group
    # high_count = (survival_tmb['tmb_group'] == 'High TMB (≥10)').sum()
    # low_count = (survival_tmb['tmb_group'] == 'Low TMB (<10)').sum()
    
    # TODO: Print group info
    # print_info(f"High TMB group: {high_count} patients")
    # print_info(f"Low TMB group: {low_count} patients")
    # print_info(f"Total patients with complete data: {len(survival_tmb)}")
    
    # I'LL PROVIDE THE KAPLAN-MEIER PLOTTING CODE
    
    pass

def create_interactive_dashboard(integrated_df):
    pass


# Main execution
if __name__ == "__main__":
    # Test loading data
    print_step(1, 6, "Loading processed data")
    tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df = load_processed_data()
    
    print_header("Data loaded successfully - ready for visualization!")
    
    print_step(2, 6, "Creating gene frequency plot")
    plot_gene_frequency(top_genes_df)
    
    print_step(3, 6, "Creating TMB distribution plot")
    plot_tmb_distribution(tmb_df)
    
    print_step(4, 6, "Creating oncoplot")
    plot_oncoplot(mutation_matrix, top_genes_df)
    
    # print_step(5, 6, "Creating survival curves")
    # plot_survival_curves(survival_df, tmb_df)
    
    # print_step(6, 6, "Creating interactive dashboard")
    # create_interactive_dashboard(integrated_df)