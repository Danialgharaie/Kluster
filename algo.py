from itertools import combinations_with_replacement as cwr
from multiprocessing import Pool
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
from umap import UMAP

from align import run_alignment


def compute_distance_matrix(
    proteins: Dict[str, str],
    alignment_tool: str,
    use_tmscore: bool,
    use_rmsd: bool,
    num_processes: int,
) -> np.ndarray:
    """Compute pairwise distance matrix using multiprocessing."""
    protein_ids = sorted(proteins.keys())
    combos = list(cwr(range(len(protein_ids)), 2))

    with Pool(processes=num_processes) as pool:
        scores = list(
            tqdm(
                pool.imap(
                    lambda c: run_alignment(
                        proteins[protein_ids[c[0]]],
                        proteins[protein_ids[c[1]]],
                        alignment_tool,
                        use_tmscore,
                        use_rmsd,
                    ),
                    combos,
                ),
                total=len(combos),
                desc="Aligning structures",
            )
        )

    # Build distance matrix
    n = len(protein_ids)
    matrix = np.zeros((n, n))
    for idx, (i, j) in enumerate(combos):
        for metric, value in scores[idx].items():
            matrix[i, j] = matrix[j, i] = value

    return matrix


def perform_clustering(
    matrix: np.ndarray,
    method: str,
    dimensions: int,
    scale: bool,
    perplexity: float = 30.0,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
) -> np.ndarray:
    """Perform dimensionality reduction and visualization."""
    # Handle missing values
    imputer = SimpleImputer(strategy="mean")
    matrix_imputed = imputer.fit_transform(matrix)

    # Scale features if requested
    matrix_processed = matrix_imputed
    if scale:
        scaler = StandardScaler()
        matrix_processed = scaler.fit_transform(matrix_imputed)

    # Perform dimensionality reduction
    if method == "TSNE":
        proj = TSNE(
            n_components=dimensions,
            perplexity=perplexity,
            random_state=42,
        ).fit_transform(matrix_processed)
    elif method == "PCA":
        proj = PCA(
            n_components=dimensions,
            random_state=42,
        ).fit_transform(matrix_processed)
    elif method == "UMAP":
        proj = UMAP(
            n_components=dimensions,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=42,
        ).fit_transform(matrix_processed)

    return proj


def visualize_projection(
    proj: np.ndarray,
    protein_ids: List[str],
    output_file: str,
    method: str,
    dimensions: int,
) -> None:
    """Generate 2D/3D visualization of the projection."""
    if dimensions == 2:
        plt.figure(figsize=(10, 10))
        plt.scatter(proj[:, 0], proj[:, 1])
        for i, txt in enumerate(protein_ids):
            plt.annotate(txt, (proj[i, 0], proj[i, 1]))
    else:  # 3D
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(proj[:, 0], proj[:, 1], proj[:, 2])
        for i, txt in enumerate(protein_ids):
            ax.text(proj[i, 0], proj[i, 1], proj[i, 2], txt)

    plt.title(f"{method} projection of protein structures")
    plt.savefig(output_file)
    plt.close()
