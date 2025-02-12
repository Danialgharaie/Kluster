import argparse
import glob
import os

import pandas as pd

from algo import (
    compute_distance_matrix,
    reduce_dimensions,
    cluster_projection,
    visualize_projection,
)
from .align import check_alignment_tool


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Cluster protein structures based on structural similarity metrics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="pdbs/",
        help="Directory containing PDB files",
    )
    parser.add_argument(
        "--alignment-tool",
        type=str,
        choices=["TMAlign", "USalign"],
        default="TMAlign",
        help="Structural alignment tool to use",
    )
    parser.add_argument(
        "--no-tmscore",
        dest="use_tmscore",
        action="store_false",
        help="Disable TM-score calculation",
    )
    parser.add_argument(
        "--no-rmsd",
        dest="use_rmsd",
        action="store_false",
        help="Disable RMSD calculation",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["UMAP", "TSNE", "PCA"],
        default="UMAP",
        help="Dimensionality reduction method",
    )
    parser.add_argument(
        "--dimensions",
        type=int,
        choices=[2, 3],
        default=2,
        help="Number of dimensions for projection",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of parallel processes for alignment",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="projection.png",
        help="Output file for the projection plot",
    )
    parser.add_argument(
        "--matrix-out",
        type=str,
        default="projection_matrix.tsv",
        help="Output file for the projection coordinates (TSV format)",
    )
    parser.add_argument(
        "--scale",
        action="store_true",
        default=False,
        help="Scale features before dimensionality reduction",
    )
    # UMAP-specific parameters
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for UMAP",
    )
    parser.add_argument(
        "--min-dist",
        type=float,
        default=0.1,
        help="Minimum distance for UMAP",
    )
    # t-SNE specific parameters
    parser.add_argument(
        "--perplexity",
        type=float,
        default=30.0,
        help="Perplexity parameter for t-SNE",
    )
    # DBSCAN parameters
    parser.add_argument(
        "--eps",
        type=float,
        default=0.5,
        help="DBSCAN eps parameter",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=5,
        help="DBSCAN min_samples parameter",
    )

    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    # Validate input
    if not any([args.use_tmscore, args.use_rmsd]):
        raise ValueError("At least one of TM-score or RMSD must be enabled")

    # Check alignment tool
    check_alignment_tool(args.alignment_tool)

    # Load PDB files
    proteins = {}
    for pdb_file in glob.glob(os.path.join(args.input_dir, "*.pdb")):
        protein_id = os.path.splitext(os.path.basename(pdb_file))[0]
        proteins[protein_id] = pdb_file

    if not proteins:
        raise ValueError(f"No PDB files found in {args.input_dir}")

    # Compute distance matrix and get sorted protein IDs
    matrix, protein_ids = compute_distance_matrix(
        proteins=proteins,
        alignment_tool=args.alignment_tool,
        use_tmscore=args.use_tmscore,
        use_rmsd=args.use_rmsd,
        num_processes=args.processes,
    )

    # Perform dimensionality reduction
    proj = reduce_dimensions(
        matrix=matrix,
        method=args.method,
        dimensions=args.dimensions,
        scale=args.scale,
        perplexity=args.perplexity,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
    )

    # Cluster the projection
    cluster_labels = cluster_projection(
        proj=proj,
        eps=args.eps,
        min_samples=args.min_samples,
    )

    # Save projection coordinates with cluster labels
    df_proj = pd.DataFrame(
        proj,
        index=protein_ids,  # Use the same sorted protein_ids
        columns=[f"component_{i+1}" for i in range(args.dimensions)]
    )
    df_proj['cluster'] = cluster_labels
    df_proj.to_csv(args.matrix_out, sep='\t')
    print(f"Projection coordinates saved to: {args.matrix_out}")

    # Generate visualization
    visualize_projection(
        proj=proj,
        protein_ids=protein_ids,
        output_file=args.output,
        method=args.method,
        dimensions=args.dimensions,
        cluster_labels=cluster_labels,
    )
    print(f"Projection plot saved to: {args.output}")


if __name__ == "__main__":
    main()
