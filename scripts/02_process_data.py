"""
Data Processing Script for TCGA LUAD Analysis
==============================================
This script takes the raw mutation and clinical data downloaded from cBioPortal
and processes it into analysis-ready formats. 

Key concepts:
- Tumor Mutation Burden (TMB): Total number of mutations per sample, normalized 
  by sequencing coverage. Higher TMB often correlates with better immunotherapy response.
- Driver vs Passenger mutations: Driver mutations contribute to cancer development,
  while passenger mutations are just along for the ride.
- Data integration: Combining mutation data with clinical outcomes for analysis.
"""

import pandas as pd
import numpy as np
import os
import sys
import json

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.display import (
    print_header, 
    print_step, 
    print_success, 
    print_info,
    print_summary,
    print_dataframe_info
)

RAW_DATA_DIR = 'data/raw'
PROCESSED_DATA_DIR = 'data/processed'
os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)

print_header("TCGA LUAD Data Processing")

def load_raw_data():
    """
    Load the three raw CSV files we downloaded.
    
    ... (rest of your docstring)
    """
    mutations_df = pd.read_csv(f'{RAW_DATA_DIR}/mutations.csv')
    clinical_sample_df = pd.read_csv(f'{RAW_DATA_DIR}/clinical_sample.csv')
    clinical_patient_df = pd.read_csv(f'{RAW_DATA_DIR}/clinical_patient.csv')

    print_dataframe_info(mutations_df, "Mutations")
    print_success(f"Loaded {len(mutations_df)} mutations")
    print('-' * 70)
    print_dataframe_info(clinical_sample_df, "Clinical Samples")
    print_success(f"Loaded {len(clinical_sample_df)} clinical samples")
    print('-' * 70)
    print_dataframe_info(clinical_patient_df, "Patients")
    print_success(f"Loaded {len(clinical_patient_df)} patients")
    print('-' * 70)

    return mutations_df, clinical_sample_df, clinical_patient_df

def extract_gene_symbols(mutations_df):
    """
    Extract gene symbols from the nested 'gene' column.
    
    Args:
        mutations_df (pd.DataFrame): Raw mutations data
        
    Returns:
        pd.DataFrame: mutations_df with new 'gene_symbol' column
    """

    # gene example: "{'entrezGeneId': 18, 'hugoGeneSymbol': 'ABAT', 'type': 'protein-coding'}"

    mutations_df['gene_symbol'] = mutations_df['gene'].apply(
        lambda x: json.loads(x.replace("'", '"'))['hugoGeneSymbol'])





# Main execution
if __name__ == "__main__":
    print_step(1, 7, "Loading raw data")
    mutations_df, clinical_sample_df, clinical_patient_df = load_raw_data()
    
    # Save processed data
    print("="*70)
    print("Extracting Gene Symbols...")
    print("="*70)

    # Extract gene symbols from mutations dataframe
    extract_gene_symbols(mutations_df)

    mutations_df.to_csv('data/processed/mutations.csv', index=False)
    print(f"✓ Saved: data/raw/mutations.csv ({len(mutations_df)} rows)")

    # print_step(2, 7, "Extracting gene symbols")
    # mutations_df = extract_gene_symbols(mutations_df)
    
    # print_step(3, 7, "Calculating tumor mutation burden")
    # tmb_df = calculate_tmb(mutations_df)
    
    # print_step(4, 7, "Identifying top mutated genes")
    # top_genes_df = identify_top_mutated_genes(mutations_df, top_n=20)
    
    # print_step(5, 7, "Creating mutation matrix")
    # mutation_matrix = create_mutation_matrix(mutations_df, top_genes_df['gene_symbol'].tolist())
    
    # print_step(6, 7, "Preparing survival data")
    # survival_df = prepare_survival_data(clinical_patient_df)
    
    # print_step(7, 7, "Integrating datasets")
    # integrated_df = integrate_datasets(tmb_df, mutation_matrix, survival_df, clinical_sample_df)
    
    # print_header("Saving processed data")
    # save_processed_data(tmb_df, mutation_matrix, survival_df, integrated_df, top_genes_df)
    
    # print_summary("Processing Complete", {
    #     "Processed files": "5 CSV files",
    #     "Location": "data/processed/",
    #     "Next step": "Run scripts/03_visualize.py"
    # })
