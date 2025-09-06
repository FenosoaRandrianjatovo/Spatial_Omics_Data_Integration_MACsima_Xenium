"""
Xenium Explorer output generation for integrated MACsima-Xenium data.

This script creates Xenium Explorer-compatible files for visualization,
including support for SpatialGlue clustering results.
"""

import pandas as pd
from pathlib import Path
import anndata as ad
import sopa

from config import *
from utils import validate_outputs


def create_combined_table(sdata_macs, clustering_results=None):
    """
    Create a combined table with genes and proteins for Xenium Explorer.
    
    Parameters
    ----------
    sdata_macs : spatialdata.SpatialData
        Processed MACsima data with aggregated tables
    clustering_results : anndata.AnnData, optional
        AnnData object containing clustering results from SpatialGlue
        
    Returns
    -------
    anndata.AnnData
        Combined table ready for Xenium Explorer
    """
    print("[INFO] Creating combined table for Xenium Explorer...")
    
    if GENES_TABLE_KEY not in sdata_macs.tables:
        raise ValueError(f"Genes table '{GENES_TABLE_KEY}' not found")
    if PROTEINS_TABLE_KEY not in sdata_macs.tables:
        raise ValueError(f"Proteins table '{PROTEINS_TABLE_KEY}' not found")
    
    genes_table = sdata_macs.tables[GENES_TABLE_KEY]
    proteins_table = sdata_macs.tables[PROTEINS_TABLE_KEY]
    
    # Start with genes table (genes will be in X matrix for visualization)
    combined_table = genes_table.copy()
    
    # Add proteins to obsm for analysis (Explorer focuses on X matrix)
    combined_table.obsm["proteins"] = proteins_table.X
    
    # Add clustering results if provided
    if clustering_results is not None:
        print("[INFO] Adding clustering results to combined table...")
        
        # Find common cells
        common_cells = combined_table.obs_names.intersection(clustering_results.obs_names)
        
        if len(common_cells) == 0:
            print("[WARNING] No common cells found between tables and clustering results")
        else:
            # Subset both objects to common cells
            combined_table = combined_table[common_cells]
            clustering_subset = clustering_results[common_cells]
            
            # Add clustering columns as categorical
            clustering_columns = []
            
            if 'SpatialGlue_Mics_Xenium' in clustering_subset.obs.columns:
                combined_table.obs['SpatialGlue_clusters'] = clustering_subset.obs['SpatialGlue_Mics_Xenium'].astype('category')
                clustering_columns.append('SpatialGlue_clusters')
                
            if 'mclust' in clustering_subset.obs.columns:
                combined_table.obs['mclust_clusters'] = clustering_subset.obs['mclust'].astype('category')
                clustering_columns.append('mclust_clusters')
                
            print(f"[INFO] Added clustering columns: {clustering_columns}")
    
    return combined_table


def write_xenium_explorer_basic(sdata_macs, output_path=None):
    """
    Write basic Xenium Explorer files without clustering results.
    
    Parameters
    ----------
    sdata_macs : spatialdata.SpatialData
        Processed MACsima data
    output_path : Path, optional
        Custom output path
    """
    if output_path is None:
        output_path = XENIUM_EXPLORER_OUTPUT
    
    print(f"[INFO] Writing basic Xenium Explorer files to: {output_path}")
    
    # Create combined table
    combined_table = create_combined_table(sdata_macs)
    
    # Add to SpatialData
    sdata_macs.tables[COMBINED_TABLE_KEY] = combined_table
    
    # Write Xenium Explorer files
    sopa.io.write.xenium_explorer(
        sdata=sdata_macs,
        output_path=output_path,
        table_name=COMBINED_TABLE_KEY,
        shapes_key=CELLS_KEY,
        gene_column="feature_name",
        cell_id_column=None,
        image_key=IMAGE_KEY
    )
    
    return output_path


def write_xenium_explorer_with_clustering(sdata_macs, clustering_results, output_path=None):
    """
    Write Xenium Explorer files with SpatialGlue clustering results.
    
    Parameters
    ----------
    sdata_macs : spatialdata.SpatialData
        Processed MACsima data
    clustering_results : anndata.AnnData
        AnnData object with clustering results
    output_path : Path, optional
        Custom output path
    """
    if output_path is None:
        output_path = XENIUM_EXPLORER_SPATIALGLUE
    
    print(f"[INFO] Writing Xenium Explorer files with clustering to: {output_path}")
    
    # Create combined table with clustering
    combined_table = create_combined_table(sdata_macs, clustering_results)
    
    # Add to SpatialData
    sdata_macs.tables[COMBINED_TABLE_KEY] = combined_table
    
    # Determine which observation columns to include
    obs_columns = []
    if 'SpatialGlue_clusters' in combined_table.obs.columns:
        obs_columns.append('SpatialGlue_clusters')
    if 'mclust_clusters' in combined_table.obs.columns:
        obs_columns.append('mclust_clusters')
    
    # Write using SOPA's explorer writer
    sopa.io.explorer.write(
        path=str(output_path),
        sdata=sdata_macs,
        table_key=COMBINED_TABLE_KEY,
        shapes_key=CELLS_KEY,
        image_key=IMAGE_KEY,
        gene_column="feature_name",
        save_h5ad=True,
        run_name="MACSima ROI1 with SpatialGlue"
    )

    return output_path


def main_basic_export(sdata_macs):
    """Export basic Xenium Explorer files."""
    print("=" * 60)
    print("XENIUM EXPLORER EXPORT - BASIC")
    print("=" * 60)
    
    try:
        output_path = write_xenium_explorer_basic(sdata_macs)
        
        if Path(output_path).exists():
            print(f"[SUCCESS] Xenium Explorer files created at: {output_path}")
        else:
            print(f"[ERROR] Failed to create Xenium Explorer files at: {output_path}")
            
        return output_path
        
    except Exception as e:
        print(f"[ERROR] Export failed: {str(e)}")
        raise


def main_clustering_export(sdata_macs, clustering_results):
    """Export Xenium Explorer files with clustering results."""
    print("=" * 60)
    print("XENIUM EXPLORER EXPORT - WITH CLUSTERING")
    print("=" * 60)
    
    try:
        output_path = write_xenium_explorer_with_clustering(sdata_macs, clustering_results)
        
        if Path(output_path).exists():
            print(f"[SUCCESS] Xenium Explorer files with clustering created at: {output_path}")
        else:
            print(f"[ERROR] Failed to create Xenium Explorer files at: {output_path}")
            
        return output_path
        
    except Exception as e:
        print(f"[ERROR] Export with clustering failed: {str(e)}")
        raise


if __name__ == "__main__":
    # This would typically be run after preprocessing and SpatialGlue analysis
    print("[INFO] This script should be run after preprocessing.py and spatialglue_integration.py")
    print("[INFO] Import the functions and call them with your processed data.")
