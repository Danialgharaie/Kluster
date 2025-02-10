# Kluster: Protein Structure Clustering Tool

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
                 --method UMAP \      # UMAP, TSNE, or PCA
                 --dimensions 2 \     # 2 or 3
                 --processes 4 \      # Number of parallel processes
                 --scale             # Scale features before reduction
```

### Command Line Arguments

- `--input-dir`: Directory containing PDB files (default: "pdbs/")
- `--alignment-tool`: Structural alignment tool to use ["TMAlign", "USalign"] (default: "TMAlign")
- `--no-tmscore`: Disable TM-score calculation
- `--no-rmsd`: Disable RMSD calculation
- `--method`: Dimensionality reduction method ["UMAP", "TSNE", "PCA"] (default: "UMAP")
- `--dimensions`: Number of dimensions for projection [2, 3] (default: 2)
- `--processes`: Number of parallel processes for alignment (default: 1)
- `--output`: Output file for the projection plot (default: "projection.png")
- `--scale`: Scale features before dimensionality reduction

Method-specific parameters:
- UMAP:
  - `--n-neighbors`: Number of neighbors (default: 15)
  - `--min-dist`: Minimum distance (default: 0.1)
- t-SNE:
  - `--perplexity`: Perplexity parameter (default: 30.0)

## Output

The tool generates a visualization plot showing the clustering of protein structures in either 2D or 3D space. The plot is saved as an image file (default: "projection.png").

## Code Organization

- `kluster.py`: Main script and argument parsing
- `align.py`: Structural alignment functions
- `algo.py`: Clustering and visualization functions

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

This clustering algorithm is adapted from:

Amani, K., Shivnauth, V., & Castroverde, C. D. M. (2023). CBP60‐DB: An AlphaFold‐predicted plant kingdom‐wide database of the CALMODULIN‐BINDING PROTEIN 60 protein family with a novel structural clustering algorithm. *Plant Direct*, 7(7). https://doi.org/10.1002/pld3.531