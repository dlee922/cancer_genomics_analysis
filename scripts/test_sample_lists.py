"""
Quick test to find sample lists
"""
import requests

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "luad_tcga_pan_can_atlas_2018"

url = f"{BASE_URL}/studies/{STUDY_ID}/sample-lists"
response = requests.get(url)

if response.status_code == 200:
    sample_lists = response.json()
    print("Available sample lists:")
    print("="*70)
    for sl in sample_lists:
        print(f"ID: {sl['sampleListId']}")
        print(f"   Name: {sl['name']}")
        print(f"   Sample count: {len(sl.get('sampleIds', []))}")
        print()
else:
    print(f"Error: {response.status_code}")