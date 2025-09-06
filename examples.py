#!/usr/bin/env python3
"""
Example script demonstrating the complete MACsima-Xenium integration workflow.

This script shows how to use the pipeline components step by step.
"""

import sys
from pathlib import Path

def run_example_workflow():
    """
    Example workflow demonstrating all pipeline components.
    """
    print(" MACsima-Xenium Integration Pipeline Example")
    print("=" * 60)
    
    # Example 1: Complete pipeline with default settings
    print("\n Example 1: Complete Pipeline (Default Settings)")
    print("-" * 50)
    print("Command: python complete_pipeline.py")
    print("This runs:")
    print("  • Data preprocessing")
    print("  • SpatialGlue integration")
    print("  • 6-cluster mclust analysis")
    print("  • Visualization generation")
    print("  • Xenium Explorer export")
    
    # Example 2: Custom clustering parameters
    print("\n Example 2: Custom Clustering Parameters")
    print("-" * 50)
    print("Command: python complete_pipeline.py --n-clusters 8 --clustering-method leiden")
    print("This runs the complete pipeline with:")
    print("  • 8 clusters instead of 6")
    print("  • Leiden clustering instead of mclust")
    
    # Example 3: Skip preprocessing
    print("\n Example 3: Skip Preprocessing (Data Exists)")
    print("-" * 50)
    print("Command: python complete_pipeline.py --skip-preprocessing")
    print("This skips preprocessing if data files already exist")
    
    # Example 4: SpatialGlue only
    print("\n Example 4: SpatialGlue Integration Only")
    print("-" * 50)
    print("Command: python complete_pipeline.py --mode spatialglue-only")
    print("This runs only SpatialGlue integration (requires preprocessed data)")
    
    # Example 5: Step by step
    print("\n Example 5: Step-by-Step Execution")
    print("-" * 50)
    print("# Step 1: Preprocessing only")
    print("python preprocessing.py")
    print()
    print("# Step 2: SpatialGlue integration")
    print("python spatialglue_integration.py")
    print()
    print("# Step 3: Xenium Explorer export")
    print("python -c \"from xenium_explorer_export import *; import sopa; import anndata as ad;\"")
    print("python -c \"sdata = sopa.io.macsima('/path/to/data'); adata = ad.read_h5ad('results.h5ad'); main_clustering_export(sdata, adata)\"")
    
    print("\n Key Parameters to Adjust:")
    print("-" * 50)
    print("• n_clusters: Number of clusters (default: 6)")
    print("• clustering_method: mclust, leiden, or louvain (default: mclust)")
    print("• random_seed: For reproducibility (default: 2025)")
    print("• Data paths in config.py")
    print("• Cellpose parameters (diameter, GPU usage)")
    
    print("\n Expected Outputs:")
    print("-" * 50)
    print("• RNA and protein AnnData files (.h5ad)")
    print("• Integrated data with SpatialGlue embeddings")
    print("• Clustering results (mclust/leiden/louvain)")
    print("• Visualization plots (UMAP, spatial, weights)")
    print("• Xenium Explorer files for interactive visualization")
    
    print("\n To get started, run:")
    print("python complete_pipeline.py")
    print("\nOr for help:")
    print("python complete_pipeline.py --help")


def check_requirements():
    """
    Check if all required packages are installed.
    """
    print("\n Checking Requirements...")
    print("-" * 30)
    
    required_packages = [
        ('numpy', 'numpy'),
        ('pandas', 'pandas'),
        ('scanpy', 'scanpy'),
        ('anndata', 'anndata'),
        ('sopa', 'sopa'),
        ('SpatialGlue', 'SpatialGlue'),
        ('torch', 'torch'),
        ('rpy2', 'rpy2'),
    ]
    
    missing_packages = []
    
    for package_name, import_name in required_packages:
        try:
            __import__(import_name)
            print(f" {package_name}")
        except ImportError:
            print(f" {package_name} - MISSING")
            missing_packages.append(package_name)
    
    if missing_packages:
        print(f"\n  Missing packages: {', '.join(missing_packages)}")
        print("Install with: pip install -r requirements.txt")
    else:
        print("\n All required packages are installed!")
    
    return len(missing_packages) == 0


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="MACsima-Xenium Pipeline Examples")
    parser.add_argument("--check-requirements", action="store_true", 
                       help="Check if all required packages are installed")
    
    args = parser.parse_args()
    
    if args.check_requirements:
        check_requirements()
    else:
        run_example_workflow()
