# Kluster: Protein Structure Clustering Tool

Kluster is a Python tool for clustering protein structures (monomers) based on their structural similarity. It uses structural alignment tools (TMAlign/USalign) to compute pairwise similarities between proteins and provides various dimensionality reduction methods for visualization.

## Features

- Multiple structural comparison metrics (TM-score, RMSD)
- Support for both TMalign and USalign
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

Place the compiled `TMalign` and `USalign` executables in the `bin/` directory before running the script. These tools are required for structural alignment.

```bash
# TMalign
wget https://zhanggroup.org/TM-align/TMalign.cpp
g++ -static -O3 -ffast-math -lm -o bin/TMalign TMalign.cpp
chmod +x bin/TMalign

# USalign
wget https://zhanggroup.org/US-align/bin/module/USalign.cpp
g++ -static -O3 -ffast-math -lm -o bin/USalign USalign.cpp
chmod +x bin/USalign
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
- `--alignment-tool`: Tool for structural alignment ("TMalign" or "USalign", default: "TMalign")
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

## Using in Neurosnap

If you don't want to set up and do the clustering and computations yourself, you can use Kluster service in Neurosnap to handle it all in an easy-to-use UI.

Service: <https://neurosnap.ai/service/Kluster>
Our blog post on protein clustering: <https://neurosnap.ai/blog/post/clustering-proteins-by-structure-the-ultimate-guide/67c238ea7748b987ae26d034>

## Citation

This clustering algorithm is adapted from:

Amani, K., Shivnauth, V., & Castroverde, C. D. M. (2023). CBP60‐DB: An AlphaFold‐predicted plant kingdom‐wide database of the CALMODULIN‐BINDING PROTEIN 60 protein family with a novel structural clustering algorithm. *Plant Direct*, 7(7). https://doi.org/10.1002/pld3.509

Alignment tools:

- Y. Zhang, J. Skolnick, TM-align: A protein structure alignment algorithm based on TM-score, Nucleic Acids Research, 33: 2302-2309 (2005) 
- Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang. US-align: Universal Structure Alignment of Proteins, Nucleic Acids and Macromolecular Complexes. Nature Methods, 19: 1109-1115 (2022)

Dimensionality reduction:

- PCA: Karl  Pearson  F.R.S. . (1901). LIII. On lines and planes of closest fit to systems of points in space. The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science, 2(11), 559–572. https://doi.org/10.1080/14786440109462720
- UMAP: McInnes et al., (2018). UMAP: Uniform Manifold Approximation and Projection. Journal of Open Source Software, 3(29), 861, https://doi.org/10.21105/joss.00861
- t-SNE: van der Maaten, L., & Hinton, G. (2008). Visualizing Data using t-SNE. Journal of Machine Learning Research, 9(86), 2579–2605. Retrieved from http://jmlr.org/papers/v9/vandermaaten08a.html

Clustering:
- DBSCAN: Ester, M., Kriegel, H.-P., Sander, J., Xu, X., & others. (1996). A density-based algorithm for discovering clusters in large spatial databases with noise. In kdd (Vol. 96, pp. 226–231).
