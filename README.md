# Cancer Genomics Analysis Pipeline

Analysis of TCGA cancer mutation data with ClinVar annotations.

## Setup

### Using Conda
\`\`\`bash
conda env create -f environment.yml
conda activate cancer-genomics
\`\`\`

### Using Docker
\`\`\`bash
docker-compose up
\`\`\`

## Project Structure
\`\`\`
cancer-genomics-analysis/
├── data/
│   ├── raw/           # Original TCGA/ClinVar data
│   └── processed/     # Cleaned, merged data
├── scripts/           # Python analysis scripts
├── notebooks/         # Jupyter notebooks
├── outputs/
│   ├── figures/       # Publication-ready plots
│   └── tables/        # Summary statistics
├── tests/             # Unit tests
└── docs/              # Documentation
\`\`\`

## Pipeline Steps
1. Data acquisition
2. Quality control
3. Mutation burden calculation
4. Driver gene identification
5. Survival analysis
6. Visualization