"""
Complete MACsima-Xenium integration pipeline with SpatialGlue.

This script runs the entire pipeline from preprocessing through SpatialGlue
integration and visualization export.
"""

import sys
import traceback
from pathlib import Path

def run_full_pipeline(n_clusters=6, clustering_method='mclust', 
                     random_seed=2025, skip_preprocessing=False):
    """
    Run the complete MACsima-Xenium integration pipeline.
    
    Parameters
    ----------
    n_clusters : int
        Number of clusters for SpatialGlue
    clustering_method : str
        Clustering method ('mclust', 'leiden', 'louvain')
    random_seed : int
        Random seed for reproducibility
    skip_preprocessing : bool
        Skip preprocessing if data already exists
        
    Returns
    -------
    tuple
        (sdata_macs, adata_integrated) with results
    """
    print("=" * 80)
    print("🧬 COMPLETE MACSIMA-XENIUM INTEGRATION PIPELINE WITH SPATIALGLUE")
    print("=" * 80)
    print(f"Configuration:")
    print(f"  - Clusters: {n_clusters}")
    print(f"  - Clustering method: {clustering_method}")
    print(f"  - Random seed: {random_seed}")
    print(f"  - Skip preprocessing: {skip_preprocessing}")
    print("=" * 80)
    
    try:
        from config import OUTPUT_BASE_DIR
        rna_path = OUTPUT_BASE_DIR / "adata_RNA_xenium.h5ad"
        adt_path = OUTPUT_BASE_DIR / "adata_ADT_mics.h5ad"
        
        # Step 1: Data Preprocessing (if needed)
        if not skip_preprocessing or not (rna_path.exists() and adt_path.exists()):
            print("\n🔄 STEP 1: Data Preprocessing")
            print("-" * 50)
            from preprocessing import main as run_preprocessing
            sdata_macs, adata_rna, adata_adt = run_preprocessing()
        else:
            print("\n⏭️  STEP 1: Skipping preprocessing (files exist)")
            print("-" * 50)
            # Load existing spatial data for later use
            print("[INFO] Loading existing spatial data structure...")
            import sopa
            from config import MACSIMA_DIR, XENIUM_DIR, CELLS_KEY, IMAGE_KEY
            sdata_macs = sopa.io.macsima(MACSIMA_DIR)
            sdata_xenium = sopa.io.xenium(XENIUM_DIR, cells_boundaries=True, nucleus_boundaries=True)
            
            # Add Xenium transcripts
            for key, pts in sdata_xenium.points.items():
                if key not in sdata_macs.points:
                    sdata_macs.points[key] = pts
        
        # Step 2: SpatialGlue Integration
        print("\n🔄 STEP 2: SpatialGlue Integration & Clustering")
        print("-" * 50)
        from spatialglue_integration import main_spatialglue_pipeline
        
        adata_integrated = main_spatialglue_pipeline(
            rna_path=rna_path,
            adt_path=adt_path,
            output_dir=OUTPUT_BASE_DIR,
            n_clusters=n_clusters,
            clustering_method=clustering_method,
            random_seed=random_seed
        )
        
        # Step 3: Create Xenium Explorer Export
        print("\n🔄 STEP 3: Xenium Explorer Export with Clustering")
        print("-" * 50)
        try:
            from xenium_explorer_export import main_clustering_export
            explorer_path = main_clustering_export(sdata_macs, adata_integrated)
            print(f"[SUCCESS] Xenium Explorer with clustering: {explorer_path}")
        except Exception as e:
            print(f"[WARNING] Xenium Explorer export with clustering failed: {e}")
            print("[INFO] Creating basic export instead...")
            try:
                from xenium_explorer_export import main_basic_export
                explorer_path = main_basic_export(sdata_macs)
                print(f"[SUCCESS] Basic Xenium Explorer: {explorer_path}")
            except Exception as e2:
                print(f"[ERROR] All Xenium Explorer exports failed: {e2}")
                explorer_path = "Export failed"
        
        # Final Summary
        print("\n✅ COMPLETE PIPELINE SUCCESSFULLY COMPLETED!")
        print("=" * 80)
        print("📊 RESULTS SUMMARY:")
        print(f"  📁 Output directory: {OUTPUT_BASE_DIR}")
        print(f"  🧬 RNA data: {rna_path}")
        print(f"  🧬 Protein data: {adt_path}")
        print(f"  🔗 Integrated data: {OUTPUT_BASE_DIR / f'data_Macsim_Xenium_spatialglue_cluster_{n_clusters}.h5ad'}")
        print(f"  📈 Plots directory: {OUTPUT_BASE_DIR / 'plots'}")
        print(f"  🌐 Xenium Explorer: {explorer_path}")
        print(f"  📊 Number of cells: {adata_integrated.n_obs}")
        print(f"  🧪 Number of genes: {adata_integrated.n_vars}")
        print(f"  🎯 Number of clusters: {n_clusters}")
        print("=" * 80)
        
        return sdata_macs, adata_integrated
        
    except Exception as e:
        print(f"\n❌ PIPELINE FAILED: {str(e)}")
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)


def run_spatialglue_only(rna_path=None, adt_path=None, **kwargs):
    """
    Run only the SpatialGlue integration part of the pipeline.
    
    Parameters
    ----------
    rna_path : str or Path, optional
        Path to RNA data
    adt_path : str or Path, optional
        Path to ADT data
    **kwargs
        Additional arguments for SpatialGlue pipeline
    """
    print("=" * 80)
    print("🔗 SPATIALGLUE INTEGRATION ONLY")
    print("=" * 80)
    
    try:
        from config import OUTPUT_BASE_DIR
        from spatialglue_integration import main_spatialglue_pipeline
        
        # Use default paths if not provided
        if rna_path is None:
            rna_path = OUTPUT_BASE_DIR / "adata_RNA_xenium.h5ad"
        if adt_path is None:
            adt_path = OUTPUT_BASE_DIR / "adata_ADT_mics.h5ad"
        
        # Check if files exist
        if not Path(rna_path).exists():
            raise FileNotFoundError(f"RNA data not found: {rna_path}")
        if not Path(adt_path).exists():
            raise FileNotFoundError(f"ADT data not found: {adt_path}")
        
        # Run SpatialGlue integration
        adata_integrated = main_spatialglue_pipeline(
            rna_path=rna_path,
            adt_path=adt_path,
            output_dir=OUTPUT_BASE_DIR,
            **kwargs
        )
        
        print(f"\n✅ SpatialGlue integration completed!")
        print(f"Results saved to: {OUTPUT_BASE_DIR}")
        
        return adata_integrated
        
    except Exception as e:
        print(f"\n❌ SpatialGlue integration failed: {str(e)}")
        raise


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Complete MACsima-Xenium Integration Pipeline with SpatialGlue",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Pipeline mode
    parser.add_argument(
        "--mode", 
        choices=["full", "spatialglue-only"],
        default="full",
        help="Pipeline execution mode"
    )
    
    # SpatialGlue parameters
    parser.add_argument(
        "--n-clusters", 
        type=int, 
        default=6,
        help="Number of clusters for SpatialGlue clustering"
    )
    
    parser.add_argument(
        "--clustering-method", 
        choices=["mclust", "leiden", "louvain"],
        default="mclust",
        help="Clustering method to use"
    )
    
    parser.add_argument(
        "--random-seed", 
        type=int, 
        default=2025,
        help="Random seed for reproducibility"
    )
    
    # Execution options
    parser.add_argument(
        "--skip-preprocessing", 
        action="store_true",
        help="Skip preprocessing if data files already exist"
    )
    
    # Data paths (for spatialglue-only mode)
    parser.add_argument(
        "--rna-path", 
        type=str,
        help="Path to RNA AnnData file (for spatialglue-only mode)"
    )
    
    parser.add_argument(
        "--adt-path", 
        type=str,
        help="Path to ADT AnnData file (for spatialglue-only mode)"
    )
    
    args = parser.parse_args()
    
    # Run pipeline based on mode
    if args.mode == "full":
        run_full_pipeline(
            n_clusters=args.n_clusters,
            clustering_method=args.clustering_method,
            random_seed=args.random_seed,
            skip_preprocessing=args.skip_preprocessing
        )
    elif args.mode == "spatialglue-only":
        run_spatialglue_only(
            rna_path=args.rna_path,
            adt_path=args.adt_path,
            n_clusters=args.n_clusters,
            clustering_method=args.clustering_method,
            random_seed=args.random_seed
        )
