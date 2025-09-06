"""
Complete pipeline runner for MACsima-Xenium data integration.

This script orchestrates the entire preprocessing pipeline and prepares
data for SpatialGlue integration.
"""

import sys
import traceback
from pathlib import Path

# Import pipeline modules
from preprocessing import main as run_preprocessing
from spatialglue_integration import main_spatialglue_pipeline
from xenium_explorer_export import main_basic_export, main_clustering_export


def run_complete_pipeline():
    """
    Run the complete preprocessing pipeline.
    
    Returns
    -------
    tuple
        (sdata_macs, adata_rna, adata_adt) for further analysis
    """
    print("=" * 80)
    print("COMPLETE MACSIMA-XENIUM INTEGRATION PIPELINE")
    print("=" * 80)
    
    try:
        # Step 1: Run preprocessing
        print("\nData Preprocessing")
        print("-" * 50)
        sdata_macs, adata_rna, adata_adt = run_preprocessing()
        
        # Step 2: Create basic Xenium Explorer export
        print("\nBasic Xenium Explorer Export")
        print("-" * 50)
        explorer_path = main_basic_export(sdata_macs)
        
        print("\nPIPELINE COMPLETED SUCCESSFULLY!")
        print("=" * 80)
        print("NEXT STEPS:")
        print("1. Run SpatialGlue integration using the generated AnnData objects")
        print("2. Use xenium_explorer_export.py to create visualization with clustering results")
        print(f"3. Open Xenium Explorer with: {explorer_path}")
        print("=" * 80)
        
        return sdata_macs, adata_rna, adata_adt
        
    except Exception as e:
        print(f"\n PIPELINE FAILED: {str(e)}")
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)


def run_complete_pipeline_with_spatialglue(n_clusters=6, clustering_method='mclust'):
    """
    Run the complete pipeline including SpatialGlue integration.
    
    Parameters
    ----------
    n_clusters : int
        Number of clusters for SpatialGlue
    clustering_method : str
        Clustering method to use
        
    Returns
    -------
    tuple
        (sdata_macs, adata_integrated) with all results
    """
    print("=" * 80)
    print("COMPLETE PIPELINE WITH SPATIALGLUE INTEGRATION")
    print("=" * 80)
    
    try:
        # Run preprocessing
        print("\n  Data Preprocessing")
        print("-" * 50)
        sdata_macs, adata_rna, adata_adt = run_preprocessing()
        
        # Run SpatialGlue integration
        print("\n SpatialGlue Integration")
        print("-" * 50)
        from config import OUTPUT_BASE_DIR
        rna_path = OUTPUT_BASE_DIR / "adata_RNA_xenium.h5ad"
        adt_path = OUTPUT_BASE_DIR / "adata_ADT_mics.h5ad"
        
        adata_integrated = main_spatialglue_pipeline(
            rna_path, adt_path, 
            n_clusters=n_clusters, 
            clustering_method=clustering_method
        )
        
        #  Create Xenium Explorer export with clustering
        print("\nXenium Explorer Export with Clustering")
        print("-" * 50)
        try:
            explorer_path = main_clustering_export(sdata_macs, adata_integrated)
            print(f"[SUCCESS] Xenium Explorer with clustering: {explorer_path}")
        except Exception as e:
            print(f"[WARNING] Xenium Explorer export failed: {e}")
            print("[INFO] Creating basic export instead...")
            explorer_path = main_basic_export(sdata_macs)
        
        print("\n COMPLETE PIPELINE WITH SPATIALGLUE COMPLETED!")
        print("=" * 80)
        print("RESULTS:")
        print(f"- Integrated data with clustering: {OUTPUT_BASE_DIR / f'data_Macsim_Xenium_spatialglue_cluster_{n_clusters}.h5ad'}")
        print(f"- Visualization plots: {OUTPUT_BASE_DIR / 'plots'}")
        print(f"- Xenium Explorer: {explorer_path}")
        print("=" * 80)
        
        return sdata_macs, adata_integrated
        
    except Exception as e:
        print(f"\n COMPLETE PIPELINE FAILED: {str(e)}")
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)


def run_with_clustering(clustering_results_path):
    """
    Run pipeline with existing clustering results.
    
    Parameters
    ----------
    clustering_results_path : str or Path
        Path to AnnData file with clustering results
    """
    import anndata as ad
    from xenium_explorer_export import main_clustering_export
    
    print("=" * 80)
    print("PIPELINE WITH EXISTING CLUSTERING RESULTS")
    print("=" * 80)
    
    try:
        # Load clustering results
        print(f"[INFO] Loading clustering results from: {clustering_results_path}")
        clustering_results = ad.read_h5ad(clustering_results_path)
        
        # Run preprocessing
        print("\n Data Preprocessing")
        print("-" * 50)
        sdata_macs, adata_rna, adata_adt = run_preprocessing()
        
        # Create Xenium Explorer export with clustering
        print("\n Xenium Explorer Export with Clustering")
        print("-" * 50)
        explorer_path = main_clustering_export(sdata_macs, clustering_results)
        
        print("\n PIPELINE WITH CLUSTERING COMPLETED!")
        print("=" * 80)
        print(f"Xenium Explorer files available at: {explorer_path}")
        print("=" * 80)
        
        return sdata_macs, adata_rna, adata_adt
        
    except Exception as e:
        print(f"\n PIPELINE FAILED: {str(e)}")
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="MACsima-Xenium Integration Pipeline")
    parser.add_argument(
        "--clustering-results", 
        type=str, 
        help="Path to AnnData file with clustering results (optional)"
    )
    parser.add_argument(
        "--spatialglue", 
        action="store_true",
        help="Run complete pipeline including SpatialGlue integration"
    )
    parser.add_argument(
        "--n-clusters", 
        type=int, 
        default=6,
        help="Number of clusters for SpatialGlue clustering (default: 6)"
    )
    parser.add_argument(
        "--clustering-method", 
        type=str, 
        default="mclust",
        choices=["mclust", "leiden", "louvain"],
        help="Clustering method to use (default: mclust)"
    )
    
    args = parser.parse_args()
    
    if args.clustering_results:
        if not Path(args.clustering_results).exists():
            print(f"[ERROR] Clustering results file not found: {args.clustering_results}")
            sys.exit(1)
        run_with_clustering(args.clustering_results)
    elif args.spatialglue:
        run_complete_pipeline_with_spatialglue(
            n_clusters=args.n_clusters,
            clustering_method=args.clustering_method
        )
    else:
        run_complete_pipeline()
