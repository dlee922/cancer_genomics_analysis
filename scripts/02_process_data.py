"""
Data Processing Script for TCGA LUAD Analysis
==============================================
This script takes the raw mutation and clinical data downloaded from cBioPortal
and processes it into analysis-ready formats.
"""

import pandas as pd
import os
import sys
import json

# Add parent directory to path to import utils
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.display import (
    print_header, 
    print_step, 
    print_success, 
    print_info,
    print_summary,
    print_error
)

# Setup paths
RAW_DATA_DIR = 'data/raw'
PROCESSED_DATA_DIR = 'data/processed'
os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)

print_header("TCGA LUAD Data Processing")


def load_raw_data():
    """
    Load the three raw CSV files we downloaded.
    
    Returns:
        tuple: (mutations_df, clinical_sample_df, clinical_patient_df)
    """
    print_info("Loading mutations data...")
    mutations_df = pd.read_csv(f'{RAW_DATA_DIR}/mutations.csv')
    
    print_info("Loading clinical sample data...")
    clinical_sample_df = pd.read_csv(f'{RAW_DATA_DIR}/clinical_sample.csv')
    
    print_info("Loading clinical patient data...")
    clinical_patient_df = pd.read_csv(f'{RAW_DATA_DIR}/clinical_patient.csv')
    
    # Print shapes: (rows, columns)
    print_info(f"Mutations: {mutations_df.shape[0]:,} rows, {mutations_df.shape[1]} columns")
    print_info(f"Clinical samples: {clinical_sample_df.shape[0]:,} rows")
    print_info(f"Clinical patients: {clinical_patient_df.shape[0]:,} rows")
    
    print_success("All data loaded successfully")
    
    return mutations_df, clinical_sample_df, clinical_patient_df


def extract_gene_symbols(mutations_df):
    """
    Extract gene symbols from the nested 'gene' column.
    
    Args:
        mutations_df (pd.DataFrame): Raw mutations data
        
    Returns:
        pd.DataFrame: mutations_df with new 'gene_symbol' column
    """
    print_info("Extracting gene symbols from nested data...")
    
    # Parse JSON string and extract the 'hugoGeneSymbol' field
    mutations_df['gene_symbol'] = mutations_df['gene'].apply(
        lambda x: json.loads(x.replace("'", '"'))['hugoGeneSymbol']
    )
    
    # Count unique genes
    unique_genes = mutations_df['gene_symbol'].nunique()
    print_success(f"Extracted {unique_genes:,} unique gene symbols")
    
    return mutations_df


def calculate_tmb(mutations_df):
    """
    Calculate Tumor Mutation Burden (TMB) for each sample.
    
    Args:
        mutations_df (pd.DataFrame): Mutations data with gene_symbol column
        
    Returns:
        pd.DataFrame: TMB data with columns ['sampleId', 'mutation_count', 'tmb']
    """
    print_info("Calculating tumor mutation burden...")
    
    # Group by sample and count mutations
    mutation_counts = mutations_df.groupby('sampleId').size()
    
    # Convert to DataFrame
    tmb_df = pd.DataFrame({
        'sampleId': mutation_counts.index,
        'mutation_count': mutation_counts.values
    })
    
    # Calculate TMB
    EXOME_SIZE_MB = 38
    tmb_df['tmb'] = tmb_df['mutation_count'] / EXOME_SIZE_MB
    
    # Print summary statistics
    print_info(f"TMB Statistics:")
    print_info(f"  Mean TMB: {tmb_df['tmb'].mean():.2f} mutations/Mb")
    print_info(f"  Median TMB: {tmb_df['tmb'].median():.2f} mutations/Mb")
    print_info(f"  Min TMB: {tmb_df['tmb'].min():.2f} mutations/Mb")
    print_info(f"  Max TMB: {tmb_df['tmb'].max():.2f} mutations/Mb")
    
    print_success(f"Calculated TMB for {len(tmb_df)} samples")
    
    return tmb_df


def identify_top_mutated_genes(mutations_df, top_n=20):
    """
    Identify the most frequently mutated genes across all samples.
    
    Args:
        mutations_df (pd.DataFrame): Mutations data
        top_n (int): Number of top genes to return (default: 20)
        
    Returns:
        pd.DataFrame: Top mutated genes with frequencies
    """
    print_info(f"Identifying top {top_n} mutated genes...")
    
    # Count how many samples each gene is mutated in
    gene_frequency = mutations_df.groupby('gene_symbol')['sampleId'].nunique()
    
    # Sort in descending order (most mutated first)
    gene_frequency = gene_frequency.sort_values(ascending=False)
    
    # Take top N genes
    top_genes = gene_frequency.head(top_n)
    
    # Calculate percentage of samples with this mutation
    total_samples = mutations_df['sampleId'].nunique()
    
    top_genes_df = pd.DataFrame({
        'gene_symbol': top_genes.index,
        'sample_count': top_genes.values,
        'percentage': (top_genes.values / total_samples * 100).round(2)
    })
    
    # Reset index so gene_symbol becomes a regular column
    top_genes_df = top_genes_df.reset_index(drop=True)
    
    # Print top 10
    print_info("Top 10 most frequently mutated genes:")
    for idx, row in top_genes_df.head(10).iterrows():
        print_info(f"  {row['gene_symbol']}: {row['sample_count']} samples ({row['percentage']}%)")
    
    print_success(f"Identified top {top_n} driver genes")
    
    return top_genes_df


def create_mutation_matrix(mutations_df, top_genes):
    """
    Create a binary mutation matrix (samples × genes).
        
    Args:
        mutations_df (pd.DataFrame): Mutations data
        top_genes (list): List of gene symbols to include
        
    Returns:
        pd.DataFrame: Binary mutation matrix
    """
    print_info(f"Creating mutation matrix for {len(top_genes)} genes...")
    
    # Filter to only include top genes
    filtered_mutations = mutations_df[mutations_df['gene_symbol'].isin(top_genes)]
    
    print_info(f"  Filtered to {len(filtered_mutations):,} mutations in top genes")
    
    # Create a binary indicator (1 = mutation present)
    filtered_mutations = filtered_mutations.copy()
    filtered_mutations['mutated'] = 1
    
    # Pivot to create matrix: samples as rows, genes as columns
    mutation_matrix = filtered_mutations.pivot_table(
        index='sampleId',
        columns='gene_symbol',
        values='mutated',
        aggfunc='max',
        fill_value=0
    )
    
    # Ensure all values are 0 or 1 (should already be, but good practice)
    mutation_matrix = mutation_matrix.astype(int)
    
    print_info(f"  Matrix shape: {mutation_matrix.shape[0]} samples × {mutation_matrix.shape[1]} genes")
    print_success("Mutation matrix created")
    
    return mutation_matrix


def prepare_survival_data(clinical_patient_df):
    """
    Prepare survival data for Kaplan-Meier analysis.
    
    Args:
        clinical_patient_df (pd.DataFrame): Patient clinical data
        
    Returns:
        pd.DataFrame: Survival data with columns ['patientId', 'time', 'event']
    """
    print_info("Preparing survival data...")

    print_info("Available columns:")
    print_info(f"  {', '.join(clinical_patient_df.columns.tolist())}")
    
    survival_cols = [col for col in clinical_patient_df.columns 
                     if 'OS' in col or 'SURVIVAL' in col or 'STATUS' in col]
    
    print_info(f"Found potential survival columns: {survival_cols}")
    
    # Find status and time columns in TCGA data
    time_col = None
    status_col = None
    
    for col in clinical_patient_df.columns:
        col_upper = col.upper()
        if 'OS_MONTHS' in col_upper or 'OVERALL_SURVIVAL_MONTHS' in col_upper:
            time_col = col
        if 'OS_STATUS' in col_upper or 'OVERALL_SURVIVAL_STATUS' in col_upper:
            status_col = col
    
    if time_col is None or status_col is None:
        print_error("Could not find survival columns!")
        print_info("Please manually inspect clinical_patient_df.columns and update this function")
        # Return empty DataFrame as fallback
        return pd.DataFrame(columns=['patientId', 'time', 'event'])
    
    print_info(f"Using time column: {time_col}")
    print_info(f"Using status column: {status_col}")
    
    # Create survival DataFrame
    survival_df = clinical_patient_df[['patientId', time_col, status_col]].copy()
    
    # Rename columns to standard names
    survival_df.columns = ['patientId', 'time', 'status']
    
    # Clean the data
    survival_df['time'] = pd.to_numeric(survival_df['time'], errors='coerce')
    
    # Convert status to binary (1 = dead, 0 = alive)
    if survival_df['status'].dtype == 'object':  # If it's a string
        # Extract first character if format is "1:DECEASED"
        survival_df['event'] = survival_df['status'].str.split(':').str[0]
        survival_df['event'] = pd.to_numeric(survival_df['event'], errors='coerce')
    else:
        survival_df['event'] = survival_df['status']
    
    # Remove rows with missing data
    survival_df = survival_df[['patientId', 'time', 'event']].dropna()
    
    # Ensure event is binary (0 or 1)
    survival_df['event'] = survival_df['event'].astype(int)
    
    # Print summary statistics
    print_info(f"Survival data summary:")
    print_info(f"  Total patients: {len(survival_df)}")
    print_info(f"  Events (deaths): {survival_df['event'].sum()}")
    print_info(f"  Censored (alive): {(survival_df['event'] == 0).sum()}")
    print_info(f"  Median follow-up time: {survival_df['time'].median():.1f} months")
    print_info(f"  Median survival time: {survival_df[survival_df['event']==1]['time'].median():.1f} months")
    
    print_success(f"Prepared survival data for {len(survival_df)} patients")
    
    return survival_df


def integrate_datasets(tmb_df, mutation_matrix, survival_df, clinical_sample_df):
    """
    Merge all datasets into one master DataFrame for analysis.
    
    Args:
        tmb_df (pd.DataFrame): TMB data
        mutation_matrix (pd.DataFrame): Binary mutation matrix
        survival_df (pd.DataFrame): Survival data
        clinical_sample_df (pd.DataFrame): Sample clinical data
        
    Returns:
        pd.DataFrame: Integrated dataset
    """
    print_info("Integrating all datasets...")
    
    integrated_df = tmb_df.copy()
    print_info(f"  Starting with {len(integrated_df)} samples from TMB data")
    
    # Add mutation matrix
    integrated_df = integrated_df.merge(
        mutation_matrix,
        left_on='sampleId',      # Column in integrated_df
        right_index=True,        # Index in mutation_matrix
        how='left'               # Keep all samples from TMB data
    )
    print_info(f"  After adding mutations: {integrated_df.shape}")
    
    # Extract patientId from sampleId
    integrated_df['patientId'] = integrated_df['sampleId'].str[:12]
    
    # Add survival data (patient-level)
    integrated_df = integrated_df.merge(
        survival_df,
        on='patientId',
        how='left'               # Keep all samples even if no survival data
    )
    print_info(f"  After adding survival: {integrated_df.shape}")
    print_info(f"  Samples with survival data: {integrated_df['time'].notna().sum()}")
    
    # Add clinical sample data
    useful_clinical_cols = ['sampleId']
    
    # Look for common clinical variables
    for col in clinical_sample_df.columns:
        col_upper = col.upper()
        if any(keyword in col_upper for keyword in ['STAGE', 'GRADE', 'AGE', 'GENDER', 'SEX', 'SMOKING']):
            useful_clinical_cols.append(col)
    
    if len(useful_clinical_cols) > 1:
        print_info(f"  Adding clinical variables: {', '.join(useful_clinical_cols[1:])}")
        integrated_df = integrated_df.merge(
            clinical_sample_df[useful_clinical_cols],
            on='sampleId',
            how='left'
        )
    else:
        print_info("  No additional clinical variables found")
    
    # Fill missing mutation values with 0 (no mutation)
    gene_columns = mutation_matrix.columns.tolist()
    integrated_df[gene_columns] = integrated_df[gene_columns].fillna(0).astype(int)
    
    print_info(f"Final integrated dataset: {integrated_df.shape[0]} samples × {integrated_df.shape[1]} features")
    print_success("Integration complete")
    
    return integrated_df


def save_processed_data(tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df):
    """
    Save all processed data to CSV files.
    
    Args:
        tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df: DataFrames to save
    """
    print_info("Saving processed datasets...")
    
    # Save TMB data
    tmb_df.to_csv(f'{PROCESSED_DATA_DIR}/tmb_data.csv', index=False)
    print_success(f"Saved: tmb_data.csv ({len(tmb_df)} samples)")
    
    # Save mutation matrix
    mutation_matrix.to_csv(f'{PROCESSED_DATA_DIR}/mutation_matrix.csv', index=True)
    print_success(f"Saved: mutation_matrix.csv ({mutation_matrix.shape[0]} samples × {mutation_matrix.shape[1]} genes)")
    
    # Save survival data
    survival_df.to_csv(f'{PROCESSED_DATA_DIR}/survival_data.csv', index=False)
    print_success(f"Saved: survival_data.csv ({len(survival_df)} patients)")
    
    # Save integrated dataset
    integrated_df.to_csv(f'{PROCESSED_DATA_DIR}/integrated_data.csv', index=False)
    print_success(f"Saved: integrated_data.csv ({len(integrated_df)} samples)")
    
    # Save top genes list
    top_genes_df.to_csv(f'{PROCESSED_DATA_DIR}/top_mutated_genes.csv', index=False)
    print_success(f"Saved: top_mutated_genes.csv ({len(top_genes_df)} genes)")
    
    print_info("\nAll processed files saved to: data/processed/")


# Main execution
if __name__ == "__main__":
    print_step(1, 7, "Loading raw data")
    mutations_df, clinical_sample_df, clinical_patient_df = load_raw_data()
    
    print_step(2, 7, "Extracting gene symbols")
    mutations_df = extract_gene_symbols(mutations_df)
    
    print_step(3, 7, "Calculating tumor mutation burden")
    tmb_df = calculate_tmb(mutations_df)
    
    print_step(4, 7, "Identifying top mutated genes")
    top_genes_df = identify_top_mutated_genes(mutations_df, top_n=20)
    
    print_step(5, 7, "Creating mutation matrix")
    mutation_matrix = create_mutation_matrix(mutations_df, top_genes_df['gene_symbol'].tolist())
    
    print_step(6, 7, "Preparing survival data")
    survival_df = prepare_survival_data(clinical_patient_df)
    
    print_step(7, 7, "Integrating datasets")
    integrated_df = integrate_datasets(tmb_df, mutation_matrix, survival_df, clinical_sample_df)
    
    print_header("Saving processed data")
    save_processed_data(tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df)
    
    print_summary("Processing Complete", {
        "TMB data": f"{len(tmb_df)} samples",
        "Mutation matrix": f"{mutation_matrix.shape[0]}×{mutation_matrix.shape[1]}",
        "Survival data": f"{len(survival_df)} patients",
        "Integrated data": f"{integrated_df.shape[0]} samples",
        "Top genes": f"{len(top_genes_df)} genes",
        "Output location": "data/processed/",
        "Next step": "Run scripts/03_visualize.py"
    })