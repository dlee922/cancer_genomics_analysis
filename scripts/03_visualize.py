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
    Creates a horizontal bar chart of top mutated genes.

    Args:
        top_genes_df (pd.DataFrame): Top mutated genes with frequencies
    """
    print_info("Creating gene frequency bar chart...")

    # Sort genes by percentage (frequency) and descending order 
    sorted_genes = top_genes_df.sort_values('percentage', ascending=False)

    # Grabbing only the top 15 genes so the bar chart is not humungous
    # change based on spefications
    plot_data = sorted_genes.head(15)
    
    # Prep data for matplotlib
    genes = plot_data['gene_symbol'].tolist()
    percentages = plot_data['percentage'].tolist()

    print_info(f"Plotting top {len(genes)} mutated genes")

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create horizontal bar chart
    # Reverse the order so highest is at top
    y_positions = np.arange(len(genes))
    bars = ax.barh(y_positions, percentages[::-1], color='steelblue', edgecolor='navy', linewidth=1.2)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(genes[::-1], fontsize=11)  # Reverse to match bars
    ax.set_xlabel('Percentage of Samples Mutated (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {len(genes)} Most Frequently Mutated Genes in LUAD', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add percentage labels on the bars
    for i, (bar, pct) in enumerate(zip(bars, percentages[::-1])):
        ax.text(pct + 1, bar.get_y() + bar.get_height()/2, 
                f'{pct:.1f}%', 
                va='center', fontsize=10, fontweight='bold')
    
    # Add grid to read easier
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
    Creates a histogram showing TMB distribution with clinical cutoffs.
    
    Args:
        tmb_df (pd.DataFrame): TMB data with 'tmb' column
    """
    print_info("Creating TMB distribution plot...")

    mean_tmb = tmb_df['tmb'].mean()
    median_tmb = tmb_df['tmb'].median()
    std_tmb = tmb_df['tmb'].std()
    
    # clinical cutoff
    TMB_CUTOFF = 10
    
    # Count high TMB samples
    high_tmb_count = (tmb_df['tmb'] >= TMB_CUTOFF).sum()
    total_samples = len(tmb_df)
    high_tmb_percentage = (high_tmb_count / total_samples) * 100

    print_info(f"Mean TMB: {mean_tmb:.2f} mut/Mb")
    print_info(f"Median TMB: {median_tmb:.2f} mut/Mb")
    print_info(f"High TMB samples (≥{TMB_CUTOFF}): {high_tmb_count} ({high_tmb_percentage:.1f}%)")    
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Create histogram
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
    Creates an oncoplot (heatmap) showing mutation patterns.
    
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
    sample_mutation_counts = plot_matrix.sum(axis=1)  
    sorted_samples = sample_mutation_counts.sort_values(ascending=False).index
    plot_matrix = plot_matrix.loc[sorted_samples] 
    
    # print info about what we're plotting
    print_info(f"Plotting {plot_matrix.shape[0]} samples × {plot_matrix.shape[1]} genes")
    print_info(f"Genes: {', '.join(top_15_genes[:5])}...")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(14, 10))
    
    # Create grid for layout: [gene bars, heatmap]
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 10], hspace=0.05)
    
    # Top subplot: gene frequency bars
    ax_genes = fig.add_subplot(gs[0])
    gene_freqs = top_genes_df.head(15)['percentage'].values
    ax_genes.bar(range(len(top_15_genes)), gene_freqs, color='coral', edgecolor='darkred')
    ax_genes.set_xlim(-0.5, len(top_15_genes) - 0.5)
    ax_genes.set_ylabel('Frequency\n(%)', fontsize=10, fontweight='bold')
    ax_genes.set_xticks([])
    ax_genes.spines['bottom'].set_visible(False)
    ax_genes.grid(axis='y', alpha=0.3)
    
    # Bottom subplot: main heatmap
    ax_heatmap = fig.add_subplot(gs[1])
    
    # Create the heatmap
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
    
    Args:
        survival_df (pd.DataFrame): Survival data (patientId, time, event)
        tmb_df (pd.DataFrame): TMB data (sampleId, mutation_count, tmb)
    """
    print_info("Creating Kaplan-Meier survival curves...")
    
    # Extract patientId from sampleId
    # TCGA format: "TCGA-05-4244-01" → patient is "TCGA-05-4244" (first 12 chars)
    tmb_with_patient = tmb_df.copy()
    tmb_with_patient['patientId'] = tmb_with_patient['sampleId'].str[:12]
    
    # Merge with survival data
    # Keep only patients with both TMB and survival data (inner join)
    survival_tmb = survival_df.merge(tmb_with_patient[['patientId', 'tmb']], on='patientId', how='inner')
    
    # Define TMB groups
    TMB_CUTOFF = 10
    survival_tmb['tmb_group'] = survival_tmb['tmb'].apply(
        lambda x: 'High TMB (≥10)' if x >= TMB_CUTOFF else 'Low TMB (<10)'
    )
    
    # Count patients in each group
    high_count = (survival_tmb['tmb_group'] == 'High TMB (≥10)').sum()
    low_count = (survival_tmb['tmb_group'] == 'Low TMB (<10)').sum()
    
    # Print group info
    print_info(f"High TMB group: {high_count} patients")
    print_info(f"Low TMB group: {low_count} patients")
    print_info(f"Total patients with complete data: {len(survival_tmb)}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Separate data by TMB group
    high_tmb_data = survival_tmb[survival_tmb['tmb_group'] == 'High TMB (≥10)']
    low_tmb_data = survival_tmb[survival_tmb['tmb_group'] == 'Low TMB (<10)']
    
    # Create Kaplan-Meier fitter objects for each group
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    # Fit the models
    # T = time, E = event (1=death, 0=censored), label = name for legend
    kmf_high.fit(durations=high_tmb_data['time'], 
                 event_observed=high_tmb_data['event'],
                 label=f'High TMB (≥{TMB_CUTOFF} mut/Mb, n={high_count})')
    
    kmf_low.fit(durations=low_tmb_data['time'],
                event_observed=low_tmb_data['event'],
                label=f'Low TMB (<{TMB_CUTOFF} mut/Mb, n={low_count})')
    
    # Plot the survival curves
    kmf_high.plot_survival_function(ax=ax, color='red', linewidth=2.5, ci_show=True)
    kmf_low.plot_survival_function(ax=ax, color='blue', linewidth=2.5, ci_show=True)
    
    # Perform log-rank test to check if difference is statistically significant
    # Compares the two survival curves
    from lifelines.statistics import logrank_test
    results = logrank_test(durations_A=high_tmb_data['time'],
                          durations_B=low_tmb_data['time'],
                          event_observed_A=high_tmb_data['event'],
                          event_observed_B=low_tmb_data['event'])
    
    p_value = results.p_value
    print_info(f"Log-rank test p-value: {p_value:.4f}")
    
    # Interpret p-value
    if p_value < 0.05:
        print_info("Survival difference is statistically significant (p < 0.05)")
        significance = "Significant difference"
    else:
        print_info("Survival difference is NOT statistically significant (p ≥ 0.05)")
        significance = "No significant difference"
    
    # Calculate median survival times
    median_high = kmf_high.median_survival_time_
    median_low = kmf_low.median_survival_time_
    print_info(f"Median survival - High TMB: {median_high:.1f} months")
    print_info(f"Median survival - Low TMB: {median_low:.1f} months")
    
    # Add labels and title
    ax.set_xlabel('Time (months)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Survival Probability', fontsize=13, fontweight='bold')
    ax.set_title('Kaplan-Meier Survival Curves: High vs Low TMB in LUAD\n' + 
                 f'Log-rank test: p = {p_value:.4f} ({significance})',
                 fontsize=14, fontweight='bold', pad=20)
    
    # Customize legend
    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    
    # Add grid
    ax.grid(alpha=0.3, linestyle='--')
    
    # Set y-axis to 0-1 (probability scale)
    ax.set_ylim([0, 1.05])
    
    # Add text box with statistics
    stats_text = (f'Median Survival:\n'
                 f'High TMB: {median_high:.1f} mo\n'
                 f'Low TMB: {median_low:.1f} mo\n'
                 f'p-value: {p_value:.4f}')
    
    ax.text(0.98, 0.02, stats_text,
            transform=ax.transAxes,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=10, family='monospace')
    
    # Clean up
    plt.tight_layout()
    
    # Save
    output_path = f'{OUTPUT_DIR}/survival_curves.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print_success(f"Saved: {output_path}")
    
    plt.close()

def create_interactive_dashboard(integrated_df):
    """
    Create an interactive Plotly scatter plot: TMB vs Survival Time.
    
    Args:
        integrated_df (pd.DataFrame): Integrated dataset with TMB + survival
    """
    print_info("Creating interactive TMB dashboard...")
    
    # Filter to samples with survival data
    plot_df = integrated_df[integrated_df['time'].notna()].copy()
    
    # Create readable status labels
    plot_df['status_label'] = plot_df['event'].map({1: 'Deceased', 0: 'Alive'})
    
    # Print summary
    print_info(f"Plotting {len(plot_df)} samples with complete data")
    deceased_count = (plot_df['event'] == 1).sum()
    alive_count = (plot_df['event'] == 0).sum()
    print_info(f"  Deceased: {deceased_count}, Alive/Censored: {alive_count}")
    
    # Create the interactive scatter plot
    fig = px.scatter(
        plot_df,
        x='tmb',
        y='time',
        color='status_label',
        color_discrete_map={'Deceased': '#e74c3c', 'Alive': '#3498db'},
        title='Interactive Dashboard: TMB vs Survival Time in LUAD',
        labels={
            'tmb': 'Tumor Mutation Burden (mutations/Mb)',
            'time': 'Survival Time (months)',
            'status_label': 'Status'
        },
        hover_data={
            'sampleId': True,
            'tmb': ':.2f',
            'time': ':.1f',
            'mutation_count': True,
            'status_label': True
        },
        size_max=10,
        opacity=0.7
    )
    
    # Add trendline to show correlation
    from scipy import stats
    
    # Calculate correlation for each group
    deceased_data = plot_df[plot_df['event'] == 1]
    alive_data = plot_df[plot_df['event'] == 0]
    
    if len(deceased_data) > 1:
        slope_d, intercept_d, r_d, p_d, _ = stats.linregress(deceased_data['tmb'], deceased_data['time'])
        print_info(f"Correlation (Deceased): r={r_d:.3f}, p={p_d:.4f}")
    
    if len(alive_data) > 1:
        slope_a, intercept_a, r_a, p_a, _ = stats.linregress(alive_data['tmb'], alive_data['time'])
        print_info(f"Correlation (Alive): r={r_a:.3f}, p={p_a:.4f}")
    
    # Add reference lines
    # Vertical line at TMB=10 (clinical cutoff)
    fig.add_vline(x=10, line_dash="dash", line_color="gray", 
                  annotation_text="Clinical Cutoff (10 mut/Mb)",
                  annotation_position="top")
    
    # Horizontal line at median survival
    median_survival = plot_df['time'].median()
    fig.add_hline(y=median_survival, line_dash="dash", line_color="gray",
                  annotation_text=f"Median Survival ({median_survival:.1f} mo)",
                  annotation_position="right")
    
    # Update layout for better appearance
    fig.update_layout(
        width=1200,
        height=700,
        font=dict(size=12),
        title_font_size=16,
        title_x=0.5,  # Center title
        hovermode='closest',
        template='plotly_white',
        legend=dict(
            title_text='Patient Status',
            orientation="v",
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="gray",
            borderwidth=1
        )
    )
    
    # Update axes
    fig.update_xaxes(
        title_font_size=14,
        gridcolor='lightgray',
        showgrid=True,
        zeroline=False
    )
    
    fig.update_yaxes(
        title_font_size=14,
        gridcolor='lightgray',
        showgrid=True,
        zeroline=False
    )
    
    # Add annotations with statistics
    stats_text = (
        f"Total Samples: {len(plot_df)}<br>"
        f"Deceased: {deceased_count} ({deceased_count/len(plot_df)*100:.1f}%)<br>"
        f"Alive: {alive_count} ({alive_count/len(plot_df)*100:.1f}%)<br>"
        f"Mean TMB: {plot_df['tmb'].mean():.2f} mut/Mb<br>"
        f"Median Survival: {median_survival:.1f} months"
    )
    
    fig.add_annotation(
        text=stats_text,
        xref="paper", yref="paper",
        x=0.02, y=0.98,
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="gray",
        borderwidth=1,
        borderpad=10,
        align="left",
        font=dict(size=11, family="monospace")
    )
    
    # Save as interactive HTML
    output_path = f'{OUTPUT_DIR}/interactive_dashboard.html'
    fig.write_html(output_path)
    print_success(f"Saved: {output_path}")
    print_info("Open this file in a web browser to interact with the plot!")
    
    # Optionally save as static PNG too
    try:
        # Requires kaleido: pip install kaleido
        static_path = f'{OUTPUT_DIR}/interactive_dashboard.png'
        fig.write_image(static_path, width=1200, height=700)
        print_success(f"Also saved static version: {static_path}")
    except Exception as e:
        print_info("(Static PNG not saved - install kaleido for this feature)")


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
    
    print_step(5, 6, "Creating survival curves")
    plot_survival_curves(survival_df, tmb_df)
    
    print_step(6, 6, "Creating interactive dashboard")
    create_interactive_dashboard(integrated_df)

    print_header("All Visualizations Complete!")
    print_success(f"All figures saved to: {OUTPUT_DIR}/")
    print_info("\nGenerated files:")
    print_info("  1. gene_frequency.png")
    print_info("  2. tmb_distribution.png")
    print_info("  3. oncoplot.png")
    print_info("  4. survival_curves.png")
    print_info("  5. interactive_dashboard.html (open in browser!)")