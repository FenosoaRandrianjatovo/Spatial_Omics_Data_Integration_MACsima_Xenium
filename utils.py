"""
Utility functions for MACsima and Xenium data integration pipeline.
"""

import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
from pathlib import Path


def pick_dapi_channel(sdata):
    """
    Automatically identify the DAPI channel from SpatialData object.
    
    Parameters
    ----------
    sdata : spatialdata.SpatialData
        SpatialData object containing image data
        
    Returns
    -------
    str
        Name of the DAPI channel
    """
    img_key = next(iter(sdata.images.keys()))
    img = sdata.images[img_key]
    names = []
    
    try:
        names = [str(x) for x in img.coords.get("c", [])]
    except Exception:
        pass
        
    if not names:
        return "DAPI"
        
    if "DAPI" in names:
        return "DAPI"
        
    for n in names:
        if "dapi" in n.lower():
            return n
            
    return "DAPI"


def create_anndata_from_table(table, data_type="RNA"):
    """
    Create AnnData object from SpatialData table.
    
    Parameters
    ----------
    table : anndata.AnnData
        Input table from SpatialData
    data_type : str
        Type of data ("RNA" or "ADT")
        
    Returns
    -------
    anndata.AnnData
        Processed AnnData object
    """
    # Convert to sparse matrix
    X_sparse = sparse.csr_matrix(table.X).astype(np.float32)
    
    # Create AnnData
    adata = ad.AnnData(
        X=X_sparse,
        obs=table.obs.copy(),
        var=table.var.copy()
    )
    
    # Add spatial coordinates if available
    if "spatial" in table.obsm:
        adata.obsm["spatial"] = np.asarray(
            np.rint(table.obsm["spatial"]), 
            dtype=np.int64
        )
    
    # Make names unique
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    
    return adata


def align_anndata_objects(adata_rna, adata_adt):
    """
    Align two AnnData objects to have the same cells.
    
    Parameters
    ----------
    adata_rna : anndata.AnnData
        RNA data
    adata_adt : anndata.AnnData
        ADT/protein data
        
    Returns
    -------
    tuple
        Aligned (adata_rna, adata_adt)
    """
    common_cells = adata_rna.obs_names.intersection(adata_adt.obs_names)
    
    adata_rna_aligned = adata_rna[common_cells].copy()
    adata_adt_aligned = adata_adt[common_cells].copy()
    
    # Copy spatial coordinates from RNA to ADT
    if "spatial" in adata_rna_aligned.obsm:
        adata_adt_aligned.obsm["spatial"] = adata_rna_aligned.obsm["spatial"].copy()
    
    return adata_rna_aligned, adata_adt_aligned


def save_anndata_objects(adata_rna, adata_adt, output_dir):
    """
    Save AnnData objects to specified directory.
    
    Parameters
    ----------
    adata_rna : anndata.AnnData
        RNA data to save
    adata_adt : anndata.AnnData
        ADT data to save
    output_dir : Path
        Output directory
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    rna_path = output_dir / "adata_RNA_xenium.h5ad"
    adt_path = output_dir / "adata_ADT_mics.h5ad"
    
    print(f"[INFO] Saving RNA data to: {rna_path}")
    adata_rna.write_h5ad(rna_path, compression="gzip")
    
    print(f"[INFO] Saving ADT data to: {adt_path}")
    adata_adt.write_h5ad(adt_path, compression="gzip")
    
    return rna_path, adt_path


def print_data_summary(adata_rna, adata_adt):
    """
    Print summary statistics for the processed data.
    
    Parameters
    ----------
    adata_rna : anndata.AnnData
        RNA data
    adata_adt : anndata.AnnData
        ADT data
    """
    print("=" * 50)
    print("DATA SUMMARY")
    print("=" * 50)
    print(f"RNA data shape: {adata_rna.shape}")
    print(f"Number of genes: {adata_rna.n_vars}")
    print(f"Number of cells: {adata_rna.n_obs}")
    print("-" * 30)
    print(f"ADT data shape: {adata_adt.shape}")
    print(f"Number of proteins: {adata_adt.n_vars}")
    print(f"Number of cells: {adata_adt.n_obs}")
    print("=" * 50)


def validate_outputs(output_paths):
    """
    Validate that output files were created successfully.
    
    Parameters
    ----------
    output_paths : list
        List of output file paths to check
        
    Returns
    -------
    bool
        True if all files exist
    """
    all_exist = True
    for path in output_paths:
        if not Path(path).exists():
            print(f"[WARNING] Output file not found: {path}")
            all_exist = False
        else:
            print(f"[SUCCESS] Created: {path}")
    
    return all_exist
