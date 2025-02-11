# Kluster: Protein Structure Clustering Tool (WIP)

Kluster is a Python tool for clustering protein structures based on their structural similarity. It uses structural alignment tools (TMAlign/USalign) to compute pairwise similarities between proteins and provides various dimensionality reduction methods for visualization.

## Features

- Multiple structural comparison metrics (TM-score, RMSD)
- Support for both TMAlign and USalign
- Dimensionality reduction using UMAP, t-SNE, or PCA
- DBSCAN clustering of reduced representations
- 2D/3D visualization with cluster coloring
- Parallel processing support

## Method

The tool uses pairwise structural comparisons to create a feature tensor, which is then reduced and clustered:

1. Compute n×n×m feature tensor using TM-score and/or RMSD (n proteins, m metrics)
2. Flatten to n×(n×m) matrix
3. Apply dimensionality reduction (UMAP/t-SNE/PCA) to get n×2 or n×3 projection
4. Cluster the projection using DBSCAN
5. Generate visualization with cluster-based coloring

UMAP is the default method as it generally creates more meaningful representations than PCA while being faster than t-SNE.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/Kluster.git
cd Kluster
```

2. Install dependencies:
```bash
pip install numpy pandas scikit-learn umap-learn matplotlib tqdm
```

3. Download and compile the alignment tools executables:
```bash
# TMAlign
wget https://zhanggroup.org/TM-align/TMalign.cpp
g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp

# USalign
wget https://zhanggroup.org/US-align/bin/module/USalign.cpp
g++ -static -O3 -ffast-math -lm -o USalign USalign.cpp
```

## Usage

### Basic Command

```bash
python kluster.py --input-dir pdbs/ --output plot.png
```

### Options

#### Input/Output
- `--input-dir`: Directory containing PDB files (default: "pdbs/")
- `--output`: Output file for projection plot (default: "projection.png")
- `--matrix-out`: Output file for projection coordinates (default: "projection_matrix.tsv")

#### Alignment Options
- `--alignment-tool`: Tool for structural alignment ("TMAlign" or "USalign", default: "TMAlign")
- `--no-tmscore`: Disable TM-score calculation
- `--no-rmsd`: Disable RMSD calculation
- `--processes`: Number of parallel processes (default: 1)

#### Dimensionality Reduction
- `--method`: Reduction method ("UMAP", "TSNE", or "PCA", default: "UMAP")
- `--dimensions`: Output dimensions (2 or 3, default: 2)
- `--scale`: Scale features before reduction

#### UMAP Parameters
- `--n-neighbors`: Number of neighbors (default: 15)
- `--min-dist`: Minimum distance (default: 0.1)

#### t-SNE Parameters
- `--perplexity`: Perplexity value (default: 30.0)

#### Clustering Parameters
- `--eps`: DBSCAN eps parameter (default: 0.5)
- `--min-samples`: DBSCAN min_samples parameter (default: 5)

### Outputs

The tool generates:
1. A projection plot with cluster-based coloring and legend
2. A TSV file containing the projection coordinates and cluster labels

## Example

```bash
python kluster.py --input-dir pdbs/ \
                 --output plot.png \
                 --matrix-out projection.tsv \
                 --method UMAP \
                 --dimensions 2 \
                 --eps 0.5 \
                 --min-samples 5
```

## Code Organization

- `kluster.py`: Main script and argument parsing
- `algo.py`: Core algorithms for matrix computation, dimensionality reduction, and clustering
- `align.py`: Interface to structural alignment tools

## Citation

This clustering algorithm is adapted from:

Amani, K., Shivnauth, V., & Castroverde, C. D. M. (2023). CBP60‐DB: An AlphaFold‐predicted plant kingdom‐wide database of the CALMODULIN‐BINDING PROTEIN 60 protein family with a novel structural clustering algorithm. *Plant Direct*, 7(7). https://doi.org/10.1002/pld3.531

Alignment tools:

- Y. Zhang, J. Skolnick, TM-align: A protein structure alignment algorithm based on TM-score, Nucleic Acids Research, 33: 2302-2309 (2005) 
- Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang. US-align: Universal Structure Alignment of Proteins, Nucleic Acids and Macromolecular Complexes. Nature Methods, 19: 1109-1115 (2022)
