"""
Data Processing Script for TCGA LUAD Analysis
==============================================
This script takes the raw mutation and clinical data downloaded from cBioPortal
and processes it into analysis-ready formats.
"""

import pandas as pd
import numpy as np
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
    print_dataframe_info, 
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
    # .replace("'", '"') converts Python dict format to valid JSON
    mutations_df['gene_symbol'] = mutations_df['gene'].apply(
        lambda x: json.loads(x.replace("'", '"'))['hugoGeneSymbol']
    )
    
    # Count unique genes
    # .nunique() = "number of unique values"
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
    # Think: "For each sample, how many mutations does it have?"
    mutation_counts = mutations_df.groupby('sampleId').size()
    
    # Convert to DataFrame with proper column names
    # mutation_counts.index = the sample IDs
    # mutation_counts.values = the counts
    tmb_df = pd.DataFrame({
        'sampleId': mutation_counts.index,
        'mutation_count': mutation_counts.values
    })
    
    # Calculate TMB by normalizing by exome size
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
    
    BIOLOGICAL CONTEXT:
    -------------------
    Driver genes (like TP53, KRAS) are mutated frequently because these mutations
    help cancer cells grow/survive. Passenger genes are mutated randomly and don't
    contribute to cancer - they're just "along for the ride."
    
    In lung adenocarcinoma, we expect to see:
    - TP53 (tumor suppressor, ~50% of tumors)
    - KRAS (growth signaling, ~30%)
    - STK11, KEAP1, EGFR (10-20% each)
    
    PANDAS EXPLANATION:
    -------------------
    .groupby('gene_symbol')['sampleId'].nunique()
    Breaking this down:
    1. Group mutations by gene
    2. For each gene, look at the 'sampleId' column
    3. Count how many UNIQUE samples (nunique = number of unique)
    
    Example:
    gene_symbol    sampleId
    TP53          SAMPLE1
    TP53          SAMPLE2
    TP53          SAMPLE1  (duplicate - same sample)
    KRAS          SAMPLE1
    
    After groupby:
    TP53    2  (appears in 2 samples)
    KRAS    1  (appears in 1 sample)
    
    Args:
        mutations_df (pd.DataFrame): Mutations data
        top_n (int): Number of top genes to return (default: 20)
        
    Returns:
        pd.DataFrame: Top mutated genes with frequencies
    """
    print_info(f"Identifying top {top_n} mutated genes...")
    
    # Count how many samples each gene is mutated in
    # Think: "For each gene, in how many different patients is it mutated?"
    gene_frequency = mutations_df.groupby('gene_symbol')['sampleId'].nunique()
    
    # Sort in descending order (most mutated first)
    # .sort_values(ascending=False) = sort high to low
    gene_frequency = gene_frequency.sort_values(ascending=False)
    
    # Take top N genes
    top_genes = gene_frequency.head(top_n)
    
    # Calculate percentage of samples with this mutation
    total_samples = mutations_df['sampleId'].nunique()
    
    # Create a nice DataFrame
    top_genes_df = pd.DataFrame({
        'gene_symbol': top_genes.index,
        'sample_count': top_genes.values,
        'percentage': (top_genes.values / total_samples * 100).round(2)
    })
    
    # Reset index so gene_symbol becomes a regular column
    # (otherwise it's stored as the "index" which can be confusing)
    top_genes_df = top_genes_df.reset_index(drop=True)
    
    # Print top 10 for inspection
    print_info("Top 10 most frequently mutated genes:")
    for idx, row in top_genes_df.head(10).iterrows():
        print_info(f"  {row['gene_symbol']}: {row['sample_count']} samples ({row['percentage']}%)")
    
    print_success(f"Identified top {top_n} driver genes")
    
    return top_genes_df


def create_mutation_matrix(mutations_df, top_genes):
    """
    Create a binary mutation matrix (samples × genes).
    
    BIOLOGICAL CONTEXT:
    -------------------
    This creates a "mutation profile" for each sample. Each row is a patient,
    each column is a gene, and each cell is 1 (mutated) or 0 (not mutated).
    
    This lets us visualize the mutation landscape and see patterns:
    - Which genes are mutated together? (co-occurrence)
    - Which genes are never mutated together? (mutual exclusivity)
    - Which samples have similar mutation profiles? (clustering)
    
    Example matrix:
                TP53  KRAS  EGFR
    SAMPLE1      1     1     0     (has TP53 and KRAS mutations)
    SAMPLE2      1     0     1     (has TP53 and EGFR mutations)
    SAMPLE3      0     1     0     (has only KRAS mutation)
    
    PANDAS EXPLANATION:
    -------------------
    .pivot_table() reshapes data from "long" to "wide" format
    
    Long format (what we have):
    sampleId    gene_symbol
    SAMPLE1     TP53
    SAMPLE1     KRAS
    SAMPLE2     TP53
    
    Wide format (what we want):
                TP53  KRAS
    SAMPLE1      1     1
    SAMPLE2      1     0
    
    Parameters:
    - index: what becomes rows (samples)
    - columns: what becomes column headers (genes)
    - values: what to put in cells (1 for mutation)
    - aggfunc: if duplicate entries, how to combine them
    
    Args:
        mutations_df (pd.DataFrame): Mutations data
        top_genes (list): List of gene symbols to include
        
    Returns:
        pd.DataFrame: Binary mutation matrix
    """
    print_info(f"Creating mutation matrix for {len(top_genes)} genes...")
    
    # Filter to only include top genes
    # .isin(list) checks if each value is in the provided list
    filtered_mutations = mutations_df[mutations_df['gene_symbol'].isin(top_genes)]
    
    print_info(f"  Filtered to {len(filtered_mutations):,} mutations in top genes")
    
    # Create a binary indicator (1 = mutation present)
    filtered_mutations = filtered_mutations.copy()
    filtered_mutations['mutated'] = 1
    
    # Pivot to create matrix: samples as rows, genes as columns
    mutation_matrix = filtered_mutations.pivot_table(
        index='sampleId',           # Rows = samples
        columns='gene_symbol',       # Columns = genes
        values='mutated',            # Cell values = 1
        aggfunc='max',               # If gene mutated multiple times in sample, still just 1
        fill_value=0                 # If no mutation, fill with 0
    )
    
    # Ensure all values are 0 or 1 (should already be, but good practice)
    mutation_matrix = mutation_matrix.astype(int)
    
    print_info(f"  Matrix shape: {mutation_matrix.shape[0]} samples × {mutation_matrix.shape[1]} genes")
    print_success("Mutation matrix created")
    
    return mutation_matrix


def prepare_survival_data(clinical_patient_df):
    """
    Prepare survival data for Kaplan-Meier analysis.
    
    BIOLOGICAL CONTEXT:
    -------------------
    Survival analysis answers: "How long do patients live after diagnosis?"
    
    We need two pieces of information:
    1. TIME: How long were they followed (in months)?
    2. EVENT: Did they die (1) or are they still alive/lost to follow-up (0)?
    
    "Censored" patients (event=0) are those who:
    - Are still alive at end of study
    - Were lost to follow-up
    - Died of unrelated causes
    
    This data lets us create Kaplan-Meier survival curves to compare:
    - Do TP53-mutated patients survive shorter?
    - Does high TMB predict better outcomes?
    
    PANDAS EXPLANATION:
    -------------------
    Working with messy real-world data requires:
    1. Finding the right columns (names vary by dataset)
    2. Cleaning values (converting "1:DECEASED" to 1)
    3. Handling missing data (dropna = drop rows with NaN)
    4. Type conversion (pd.to_numeric makes strings into numbers)
    
    Args:
        clinical_patient_df (pd.DataFrame): Patient clinical data
        
    Returns:
        pd.DataFrame: Survival data with columns ['patientId', 'time', 'event']
    """
    print_info("Preparing survival data...")
    
    # First, let's see what columns we have
    print_info("Available columns:")
    print_info(f"  {', '.join(clinical_patient_df.columns.tolist())}")
    
    # Look for survival-related columns
    # Common patterns: OS_MONTHS, OS_STATUS, OVERALL_SURVIVAL, DFS (disease-free survival)
    survival_cols = [col for col in clinical_patient_df.columns 
                     if 'OS' in col or 'SURVIVAL' in col or 'STATUS' in col]
    
    print_info(f"Found potential survival columns: {survival_cols}")
    
    # Try to identify time and status columns
    # This is dataset-specific - TCGA typically uses these names
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
    # Convert time to numeric (handles strings like "45.2" or "NA")
    survival_df['time'] = pd.to_numeric(survival_df['time'], errors='coerce')
    
    # Convert status to binary (1 = dead, 0 = alive)
    # TCGA format is often "1:DECEASED" or "0:LIVING"
    # We need to extract just the number
    if survival_df['status'].dtype == 'object':  # If it's a string
        # Extract first character if format is "1:DECEASED"
        survival_df['event'] = survival_df['status'].str.split(':').str[0]
        survival_df['event'] = pd.to_numeric(survival_df['event'], errors='coerce')
    else:
        survival_df['event'] = survival_df['status']
    
    # Remove rows with missing data
    # .dropna() removes any row where time or event is NaN
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
    
    BIOLOGICAL CONTEXT:
    -------------------
    Real cancer research requires integrating multiple data types:
    - Molecular (mutations, TMB)
    - Clinical (survival, stage, demographics)
    - Outcomes (did patient respond to treatment?)
    
    This integration enables questions like:
    - "In stage III patients, does TP53 mutation predict worse survival?"
    - "Is high TMB associated with better immunotherapy response?"
    
    PANDAS EXPLANATION:
    -------------------
    pd.merge() combines DataFrames based on a common column (like JOIN in SQL)
    
    Types of merges:
    - 'inner': Keep only rows that match in both DataFrames
    - 'left': Keep all rows from left DataFrame
    - 'outer': Keep all rows from both DataFrames
    
    Example:
    DataFrame A:          DataFrame B:
    sampleId  tmb        sampleId  stage
    SAMPLE1   5.2        SAMPLE1   III
    SAMPLE2   3.1        SAMPLE3   II
    
    After merge (inner):
    sampleId  tmb   stage
    SAMPLE1   5.2   III
    (SAMPLE2 and SAMPLE3 dropped - no match)
    
    Args:
        tmb_df (pd.DataFrame): TMB data
        mutation_matrix (pd.DataFrame): Binary mutation matrix
        survival_df (pd.DataFrame): Survival data
        clinical_sample_df (pd.DataFrame): Sample clinical data
        
    Returns:
        pd.DataFrame: Integrated dataset
    """
    print_info("Integrating all datasets...")
    
    # Start with TMB data
    integrated_df = tmb_df.copy()
    print_info(f"  Starting with {len(integrated_df)} samples from TMB data")
    
    # Add mutation matrix
    # mutation_matrix has sampleId as index, so we merge on index
    integrated_df = integrated_df.merge(
        mutation_matrix,
        left_on='sampleId',      # Column in integrated_df
        right_index=True,        # Index in mutation_matrix
        how='left'               # Keep all samples from TMB data
    )
    print_info(f"  After adding mutations: {integrated_df.shape}")
    
    # Extract patientId from sampleId
    # TCGA format: sampleId often contains patientId as prefix
    # Example: "TCGA-05-4244-01" → patient is "TCGA-05-4244"
    # We'll take first 12 characters (standard TCGA patient ID length)
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
    # First, let's see what useful columns are available
    useful_clinical_cols = ['sampleId']
    
    # Look for common clinical variables
    for col in clinical_sample_df.columns:
        col_upper = col.upper()
        if any(keyword in col_upper for keyword in ['STAGE', 'GRADE', 'AGE', 'GENDER', 'SEX', 'SMOKING']):
            useful_clinical_cols.append(col)
    
    if len(useful_clinical_cols) > 1:  # If we found useful columns beyond sampleId
        print_info(f"  Adding clinical variables: {', '.join(useful_clinical_cols[1:])}")
        integrated_df = integrated_df.merge(
            clinical_sample_df[useful_clinical_cols],
            on='sampleId',
            how='left'
        )
    else:
        print_info("  No additional clinical variables found")
    
    # Fill missing mutation values with 0 (no mutation)
    # Get column names that are gene symbols (from mutation matrix)
    gene_columns = mutation_matrix.columns.tolist()
    integrated_df[gene_columns] = integrated_df[gene_columns].fillna(0).astype(int)
    
    print_info(f"Final integrated dataset: {integrated_df.shape[0]} samples × {integrated_df.shape[1]} features")
    print_success("Integration complete")
    
    return integrated_df


def save_processed_data(tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df):
    """
    Save all processed data to CSV files.
    
    BIOLOGICAL CONTEXT:
    -------------------
    Saving processed data follows best practices in computational biology:
    - Raw data stays raw (never modified)
    - Processed data is saved separately
    - Reproducibility: anyone can load processed data without rerunning pipeline
    
    PANDAS EXPLANATION:
    -------------------
    .to_csv() writes a DataFrame to a CSV file
    - index=False: don't save the row numbers as a column
    - index=True: save the index (useful when index contains meaningful data)
    
    Args:
        tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df: DataFrames to save
    """
    print_info("Saving processed datasets...")
    
    # Save TMB data
    tmb_df.to_csv(f'{PROCESSED_DATA_DIR}/tmb_data.csv', index=False)
    print_success(f"Saved: tmb_data.csv ({len(tmb_df)} samples)")
    
    # Save mutation matrix (keep index since it's the sampleId)
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