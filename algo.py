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


def _run_alignment_pair(args):
    """Helper function for parallel alignment."""
    proteins, protein_ids, combo, alignment_tool, use_tmscore, use_rmsd = args
    return run_alignment(
        proteins[protein_ids[combo[0]]],
        proteins[protein_ids[combo[1]]],
        alignment_tool,
        use_tmscore,
        use_rmsd,
    )


def compute_distance_matrix(
    proteins: Dict[str, str],
    alignment_tool: str,
    use_tmscore: bool,
    use_rmsd: bool,
    num_processes: int,
) -> np.ndarray:
    """Compute pairwise distance matrix using multiprocessing.
    
    Returns:
        np.ndarray: A flattened feature matrix of shape (n, n*m) where:
            n is the number of proteins
            m is the number of features (TM-score and/or RMSD)
    """
    protein_ids = sorted(proteins.keys())
    combos = list(cwr(range(len(protein_ids)), 2))

    # Prepare arguments for parallel processing
    parallel_args = [
        (proteins, protein_ids, combo, alignment_tool, use_tmscore, use_rmsd)
        for combo in combos
    ]

    with Pool(processes=num_processes) as pool:
        scores = list(
            tqdm(
                pool.imap(_run_alignment_pair, parallel_args),
                total=len(combos),
                desc="Aligning structures",
            )
        )

    # Build feature tensors
    n = len(protein_ids)
    features = []
    
    if use_tmscore:
        tm_matrix = np.zeros((n, n))
        for idx, (i, j) in enumerate(combos):
            if 'tma_score' in scores[idx]:
                tm_matrix[i, j] = tm_matrix[j, i] = scores[idx]['tma_score']
        features.append(tm_matrix)
    
    if use_rmsd:
        rmsd_matrix = np.zeros((n, n))
        for idx, (i, j) in enumerate(combos):
            if 'rmsd' in scores[idx]:
                rmsd_matrix[i, j] = rmsd_matrix[j, i] = scores[idx]['rmsd']
        features.append(rmsd_matrix)
    
    # Stack features into a tensor and flatten
    feature_tensor = np.stack(features, axis=-1)  # Shape: (n, n, m)
    flattened_matrix = feature_tensor.reshape(n, -1)  # Shape: (n, n*m)
    
    return flattened_matrix


def perform_clustering(
    matrix: np.ndarray,
    method: str,
    dimensions: int,
    scale: bool,
    perplexity: float = 30.0,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
) -> np.ndarray:
    """Perform dimensionality reduction on the flattened feature matrix.
    
    Args:
        matrix: Flattened feature matrix of shape (n, n*m)
        method: Dimensionality reduction method (UMAP, TSNE, or PCA)
        dimensions: Output dimensions (2 or 3)
        scale: Whether to scale features before reduction
        perplexity: t-SNE perplexity parameter
        n_neighbors: UMAP n_neighbors parameter
        min_dist: UMAP min_dist parameter
    
    Returns:
        np.ndarray: Reduced dimensional representation of shape (n, dimensions)
    """
    import warnings
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=UserWarning)

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
            random_state=None,  # Allow parallelism
            n_jobs=-1,  # Use all available cores
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
    # Create figure with appropriate size
    plt.figure(figsize=(12, 8))
    
    if dimensions == 2:
        plt.scatter(proj[:, 0], proj[:, 1], alpha=0.7)
        plt.xlabel('Component 1')
        plt.ylabel('Component 2')
    else:  # 3D plot
        ax = plt.axes(projection='3d')
        ax.scatter(proj[:, 0], proj[:, 1], proj[:, 2], alpha=0.7)
        ax.set_xlabel('Component 1')
        ax.set_ylabel('Component 2')
        ax.set_zlabel('Component 3')
    
    plt.title(f"{method} projection of protein structures")
    
    # Adjust layout
    plt.tight_layout()
    
    # Save with high DPI for better quality
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
