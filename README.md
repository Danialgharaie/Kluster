# Kluster: Protein Structure Clustering Tool (WIP)

Kluster is a Python tool for clustering protein structures based on their structural similarity. It uses structural alignment tools (TMAlign/USalign) to compute pairwise similarities between proteins and provides various dimensionality reduction methods for visualization.

## Features

- Structural alignment using TMAlign or USalign
- Multiple similarity metrics (TM-score and RMSD)
- Various dimensionality reduction methods:
  - UMAP (Uniform Manifold Approximation and Projection)
  - t-SNE (t-Distributed Stochastic Neighbor Embedding)
  - PCA (Principal Component Analysis)
- Parallel processing for faster computation
- 2D and 3D visualization options
- Exports distance matrices and clustering results

## Requirements

- Python 3.6+
- TMAlign or USalign executable in the current directory
- Required Python packages:
  ```
  numpy
  pandas
  scikit-learn
  umap-learn
  matplotlib
  tqdm
  ```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/Kluster.git
   cd Kluster
   ```

2. Install required Python packages:
   ```bash
   pip install numpy pandas scikit-learn umap-learn matplotlib tqdm
   ```

3. Install the alignment tool:
   - TMAlign: Download from [Zhang Lab](https://zhanggroup.org/TM-align/)
   - USalign: Download from [Zhang Lab](https://zhanggroup.org/US-align/)
   
   Place the executable in your working directory and make it executable:
   ```bash
   chmod +x TMalign  # or USalign
   ```

## Usage

Basic usage:
```bash
python kluster.py --input-dir path/to/pdb/files
```
Advanced options:
```bash
python kluster.py --input-dir pdbs/ \
                 --alignment-tool TMAlign \
                 --method UMAP \
                 --dimensions 2 \
                 --processes 4 \
                 --output projection.png \
                 --matrix-out distances.tsv \
                 --clusters-out clusters.tsv
```

### Command Line Arguments

#### Input/Output
- `--input-dir`: Directory containing PDB files (default: "pdbs/")
- `--output`: Output file for projection plot (default: "projection.png")
- `--matrix-out`: Output file for distance matrix in TSV format (default: "distance_matrix.tsv")

#### Alignment Options
- `--alignment-tool`: Tool for structural alignment ("TMAlign" or "USalign", default: "TMAlign")
- `--no-tmscore`: Disable TM-score calculation
- `--no-rmsd`: Disable RMSD calculation
- `--processes`: Number of parallel processes for alignment (default: 1)

#### Dimensionality Reduction
- `--method`: Method for dimensionality reduction ("UMAP", "TSNE", or "PCA", default: "UMAP")
- `--dimensions`: Number of dimensions for projection (2 or 3, default: 2)
- `--scale`: Scale features before dimensionality reduction

#### Method-Specific Parameters
UMAP:
- `--n-neighbors`: Number of neighbors (default: 15)
- `--min-dist`: Minimum distance (default: 0.1)

t-SNE:
- `--perplexity`: Perplexity parameter (default: 30.0)

### Output Files

The tool generates two output files:
1. **Processed Feature Matrix** (TSV): The final processed feature matrix after imputation and scaling, ready for further analysis
2. **Visualization Plot** (PNG): 2D/3D projection of protein structures using the selected dimensionality reduction method

## Method

The tool uses a multi-step process to analyze protein structural similarity:

1. **Feature Extraction**: 
   - Computes pairwise structural comparisons using TM-Align and/or RMSD
   - Creates an n × n × m feature tensor, where:
     - n is the number of proteins
     - m is the number of features (TM-score and/or RMSD)

2. **Data Processing**:
   - Flattens the feature tensor to n × (n×m) matrix
   - Handles missing values through mean imputation
   - Optionally scales features

3. **Dimensionality Reduction**:
   - UMAP (default): Creates meaningful representations, generally outperforming PCA
   - t-SNE: Produces comparable projections but slower than UMAP
   - PCA: Available as a simpler alternative

4. **Visualization**:
   - Generates 2D/3D scatter plots
   - Saves both the processed feature matrix and visualization

## Example

```bash
python kluster.py --input-dir pdbs/ \
                 --output plot.png \
                 --matrix-out processed.tsv \
                 --method UMAP \
                 --dimensions 2
```

## Code Organization

- `kluster.py`: Main script and argument parsing
- `align.py`: Structural alignment functions
- `algo.py`: Clustering and visualization functions

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

This clustering algorithm is adapted from:

Amani, K., Shivnauth, V., & Castroverde, C. D. M. (2023). CBP60‐DB: An AlphaFold‐predicted plant kingdom‐wide database of the CALMODULIN‐BINDING PROTEIN 60 protein family with a novel structural clustering algorithm. *Plant Direct*, 7(7). https://doi.org/10.1002/pld3.531
