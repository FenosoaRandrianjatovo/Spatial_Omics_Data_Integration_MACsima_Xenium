"""
Data preprocessing pipeline for MACsima and Xenium integration.

This script handles:
1. Loading MACsima and Xenium data
2. Cell segmentation using Cellpose
3. Transcript and protein aggregation
4. Creation of separate RNA and ADT AnnData objects
"""

import os
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import spatialdata as sd
from spatialdata_io.experimental import to_legacy_anndata
import sopa

from config import *
from utils import (
    pick_dapi_channel, 
    create_anndata_from_table, 
    align_anndata_objects, 
    save_anndata_objects, 
    print_data_summary, 
    validate_outputs
)


def load_spatial_data():
    """Load MACsima and Xenium spatial data."""
    print("[INFO] Loading MACsima data...")
    sdata_macs = sopa.io.macsima(MACSIMA_DIR)
    
    print("[INFO] Loading Xenium data...")
    sdata_xenium = sopa.io.xenium(
        XENIUM_DIR, 
        cells_boundaries=True, 
        nucleus_boundaries=True
    )
    
    return sdata_macs, sdata_xenium


def setup_image_patches(sdata_macs):
    """Create image patches for processing."""
    print("[INFO] Creating image patches...")
    sopa.make_image_patches(
        sdata_macs, 
        patch_width=PATCH_WIDTH, 
        patch_overlap=PATCH_OVERLAP
    )


def perform_cell_segmentation(sdata_macs):
    """Perform cell segmentation using Cellpose if not already done."""
    print("[INFO] Checking for existing segmentation...")
    
    if CELLS_KEY not in sdata_macs.shapes:
        dapi = pick_dapi_channel(sdata_macs)
        print(f"[INFO] Using nuclear channel for Cellpose: {dapi}")
        
        print("[INFO] Running Cellpose segmentation...")
        sopa.segmentation.cellpose(
            sdata=sdata_macs,
            pretrained_model=PRETRAINED_MODEL,
            channels=dapi,
            diameter=DIAMETER_PX,
            model_type="cyto3",
            gpu=RUN_WITH_GPU,
            image_key=IMAGE_KEY,
            key_added=CELLS_KEY
        )
        print("[INFO] Cell segmentation completed!")
    else:
        print("[INFO] Using existing segmentation.")


def integrate_xenium_transcripts(sdata_macs, sdata_xenium):
    """Add Xenium transcripts to MACsima data."""
    print("[INFO] Integrating Xenium transcripts...")
    
    for key, pts in sdata_xenium.points.items():
        if key not in sdata_macs.points:
            sdata_macs.points[key] = pts
    
    print(f"[INFO] Added {len(sdata_xenium.points)} transcript datasets")


def aggregate_data(sdata_macs):
    """Aggregate both genes and proteins."""
    print("[INFO] Aggregating transcripts (genes)...")
    sopa.aggregate(
        sdata=sdata_macs,
        aggregate_genes=True,
        aggregate_channels=False,
        shapes_key=CELLS_KEY,
        key_added=GENES_TABLE_KEY,
    )

    print("[INFO] Aggregating channel intensities (proteins)...")
    sopa.aggregate(
        sdata=sdata_macs,
        aggregate_genes=False,
        aggregate_channels=True,
        image_key=IMAGE_KEY,
        shapes_key=CELLS_KEY,
        key_added=PROTEINS_TABLE_KEY,
    )


def create_anndata_objects(sdata_macs):
    """Create separate AnnData objects for RNA and ADT data."""
    print("[INFO] Building AnnData objects...")
    
    # Debug info
    print(f"[DEBUG] Available tables: {list(sdata_macs.tables.keys())}")
    
    # Build RNA AnnData from genes table
    print("[INFO] Creating RNA AnnData...")
    adata_rna = to_legacy_anndata(
        sdata_macs, 
        table_name=GENES_TABLE_KEY, 
        include_images=False
    )
    adata_rna.X = sparse.csr_matrix(adata_rna.X).astype(np.float32)
    
    if "spatial" in adata_rna.obsm:
        adata_rna.obsm["spatial"] = np.asarray(
            np.rint(adata_rna.obsm["spatial"]), 
            dtype=np.int64
        )
    
    adata_rna.obs_names_make_unique()
    adata_rna.var_names_make_unique()
    
    # Build ADT AnnData from proteins table
    print("[INFO] Creating ADT AnnData...")
    
    if PROTEINS_TABLE_KEY not in sdata_macs.tables:
        raise ValueError(f"Protein table '{PROTEINS_TABLE_KEY}' not found")
    
    protein_table = sdata_macs.tables[PROTEINS_TABLE_KEY]
    X_adt = protein_table.X
    protein_names = protein_table.var_names
    cell_names = protein_table.obs_names
    
    adata_adt = ad.AnnData(
        X=sparse.csr_matrix(X_adt).astype(np.float32),
        obs=pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=protein_names)
    )
    
    # Align objects and add spatial coordinates
    adata_rna, adata_adt = align_anndata_objects(adata_rna, adata_adt)
    
    return adata_rna, adata_adt


def main():
    """Main preprocessing pipeline."""
    print("=" * 60)
    print("MACSIMA-XENIUM DATA INTEGRATION PREPROCESSING PIPELINE")
    print("=" * 60)
    
    try:
        # Step 1: Load data
        sdata_macs, sdata_xenium = load_spatial_data()
        
        # Step 2: Setup image processing
        setup_image_patches(sdata_macs)
        print("=" * 40)
        
        # Step 3: Cell segmentation
        perform_cell_segmentation(sdata_macs)
        print("=" * 40)
        
        # Step 4: Integrate Xenium transcripts
        integrate_xenium_transcripts(sdata_macs, sdata_xenium)
        print("=" * 40)
        
        # Step 5: Aggregate data
        aggregate_data(sdata_macs)
        print("=" * 40)
        
        # Step 6: Create AnnData objects
        adata_rna, adata_adt = create_anndata_objects(sdata_macs)
        
        # Step 7: Print summary
        print_data_summary(adata_rna, adata_adt)
        
        # Step 8: Save results
        print("[INFO] Saving AnnData objects...")
        rna_path, adt_path = save_anndata_objects(adata_rna, adata_adt, OUTPUT_BASE_DIR)
        
        # Step 9: Validate outputs
        if validate_outputs([rna_path, adt_path]):
            print("[SUCCESS] Preprocessing pipeline completed successfully!")
        else:
            print("[ERROR] Some output files are missing!")
            
        # Store processed data for next steps
        return sdata_macs, adata_rna, adata_adt
        
    except Exception as e:
        print(f"[ERROR] Pipeline failed: {str(e)}")
        raise


if __name__ == "__main__":
    sdata_macs, adata_rna, adata_adt = main()
