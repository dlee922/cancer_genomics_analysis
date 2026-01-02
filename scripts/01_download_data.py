"""
Download TCGA LUAD data from cBioPortal
Uses their public API to get mutation and clinical data
"""
import requests
import json
import pandas as pd
import os
from time import sleep
import yaml

from utils.display import print_header

# create data directories if they do not exist
os.makedirs('data/raw', exist_ok=True)
os.makedirs('data/processed', exist_ok=True)

print("="*70)
print("TCGA LUAD Data Download from cBioPortal")
print("="*70)

# cBioPortal API base URL
BASE_URL = "https://www.cbioportal.org/api"

# Load configuration
with open('config.yaml', 'r') as f:
    config = yaml.safe_load(f)

STUDY_ID = config['study']['study_id']
MOLECULAR_PROFILE_ID = f"{STUDY_ID}_mutations"
SAMPLE_LIST_ID = f"{STUDY_ID}_sequenced"

print_header(f"TCGA {config['study']['cancer_abbreviation']} Data Download from cBioPortal")

def get_mutations_by_sample_list(MOLECULAR_PROFILE_ID, sample_list_id):
    """Download mutations using sample list ID"""
    print(f"\n[1/3] Downloading mutations for sample list: {sample_list_id}...")
    
    url = f"{BASE_URL}/molecular-profiles/{MOLECULAR_PROFILE_ID}/mutations"
    
    params = {
        "sampleListId": sample_list_id,
        "projection": "DETAILED"
    }
    
    response = requests.get(url, params=params)
    
    if response.status_code != 200:
        print(f"   Error: {response.status_code}")
        print(f"   Response: {response.text[:500]}")
        return pd.DataFrame()
    
    data = response.json()
    df = pd.DataFrame(data)
    print(f"   ✓ Downloaded {len(df)} mutations")
    
    return df

def get_clinical_data(study_id):
    """Download clinical data for a study"""
    print(f"\n[2/3] Downloading clinical data...")
    
    # Get all samples in the study
    url = f"{BASE_URL}/studies/{study_id}/samples"
    params = {"projection": "DETAILED"}
    
    response = requests.get(url, params=params)
    samples = response.json()
    print(f"   Found {len(samples)} samples")
    
    # Get clinical data for all samples
    url = f"{BASE_URL}/studies/{study_id}/clinical-data"
    params = {
        "clinicalDataType": "SAMPLE",
        "projection": "DETAILED"
    }
    
    response = requests.get(url, params=params)
    clinical = response.json()
    
    # Convert to DataFrame
    df = pd.DataFrame(clinical)
    print(f"   ✓ Clinical attributes: {len(df)}")
    
    # Pivot to wide format
    df_wide = df.pivot_table(
        index='sampleId',
        columns='clinicalAttributeId',
        values='value',
        aggfunc='first'
    ).reset_index()
    
    return df_wide

def get_patient_survival(study_id):
    """Download patient survival data"""
    print(f"\n[3/3] Downloading survival data...")
    
    url = f"{BASE_URL}/studies/{study_id}/patients"
    params = {"projection": "DETAILED"}
    
    response = requests.get(url, params=params)
    patients = response.json()
    
    # Get patient clinical data (includes survival)
    url = f"{BASE_URL}/studies/{study_id}/clinical-data"
    params = {
        "clinicalDataType": "PATIENT",
        "projection": "DETAILED"
    }
    
    response = requests.get(url, params=params)
    clinical = response.json()
    
    df = pd.DataFrame(clinical)
    
    # Pivot to wide format
    df_wide = df.pivot_table(
        index='patientId',
        columns='clinicalAttributeId',
        values='value',
        aggfunc='first'
    ).reset_index()
    
    print(f"   ✓ Patient records: {len(df_wide)}")
    
    return df_wide

# Main execution
if __name__ == "__main__":
    try:
        # Download data using sample list
        mutations_df = get_mutations_by_sample_list(MOLECULAR_PROFILE_ID, SAMPLE_LIST_ID)
        
        if mutations_df.empty:
            print("\n✗ No mutations downloaded. Check the error messages above.")
            exit(1)
        
        # Download clinical data
        clinical_df = get_clinical_data(STUDY_ID)
        survival_df = get_patient_survival(STUDY_ID)
        
        # Save raw data
        print("\n" + "="*70)
        print("Saving data to files...")
        print("="*70)
        
        mutations_df.to_csv('data/raw/mutations.csv', index=False)
        print(f"✓ Saved: data/raw/mutations.csv ({len(mutations_df)} rows)")
        
        clinical_df.to_csv('data/raw/clinical_sample.csv', index=False)
        print(f"✓ Saved: data/raw/clinical_sample.csv ({len(clinical_df)} rows)")
        
        survival_df.to_csv('data/raw/clinical_patient.csv', index=False)
        print(f"✓ Saved: data/raw/clinical_patient.csv ({len(survival_df)} rows)")
        
                # Print summary
        print("\n" + "="*70)
        print("Data Summary")
        print("="*70)
        print(f"Mutations: {len(mutations_df):,} variants")
        print(f"Samples: {len(clinical_df):,} samples")
        print(f"Patients: {len(survival_df):,} patients")
        
        if len(mutations_df) > 0:
            print(f"\nMutation data columns: {list(mutations_df.columns)}")
            
            # The 'gene' column contains dict objects, extract gene symbols
            if 'gene' in mutations_df.columns:
                # Extract hugoGeneSymbol from the gene dict
                if isinstance(mutations_df['gene'].iloc[0], dict):
                    print("\n   Note: 'gene' column contains nested data, extracting symbols...")
                    gene_symbols = mutations_df['gene'].apply(lambda x: x.get('hugoGeneSymbol') if isinstance(x, dict) else x)
                    print(f"Unique genes mutated: {gene_symbols.nunique():,}")
                else:
                    print(f"Unique genes mutated: {mutations_df['gene'].nunique():,}")
            
            if 'sampleId' in mutations_df.columns:
                print(f"Unique samples with mutations: {mutations_df['sampleId'].nunique():,}")
        
        print("\n✓ Data download complete!")
        print("   Files saved in data/raw/")
        print("\nNext step: Run scripts/02_process_data.py to clean and process the data")
        
    except Exception as e:
        print(f"\n✗ Error occurred: {e}")
        import traceback
        traceback.print_exc()


