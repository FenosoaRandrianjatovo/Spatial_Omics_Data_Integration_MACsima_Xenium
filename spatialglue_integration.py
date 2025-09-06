"""
SpatialGlue integration module for MACsima-Xenium data integration.

This module performs cross-modal integration using SpatialGlue and includes
clustering analysis with visualization capabilities.
"""

import os
import torch
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scipy

# SpatialGlue imports
import SpatialGlue
from SpatialGlue.preprocess import construct_neighbor_graph, clr_normalize_each_cell, pca, fix_seed
from SpatialGlue.SpatialGlue_pyG import Train_SpatialGlue
from SpatialGlue.utils import clustering

# R integration for mclust
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter

from config import OUTPUT_BASE_DIR


def setup_environment():
    """Setup R environment for mclust clustering."""
    print("[INFO] Setting up R environment for mclust...")
    
    # Clear any existing R environment variables
    for key in list(os.environ.keys()):
        if key.startswith('R_'):
            del os.environ[key]
    
    # Try to detect R installation automatically
    r_home_candidates = [
        '/usr/lib/R',
        '/usr/local/lib/R',
        '/opt/R/lib/R',
        '/Library/Frameworks/R.framework/Resources',  # macOS
        '/System/Library/Frameworks/R.framework/Resources'  # macOS system
    ]
    
    r_home = None
    for candidate in r_home_candidates:
        if os.path.exists(candidate):
            r_home = candidate
            break
    
    if r_home is None:
        print("[WARNING] Could not auto-detect R installation. Please set R_HOME manually.")
        # Fallback - try to use R from PATH
        try:
            import subprocess
            result = subprocess.run(['R', 'RHOME'], capture_output=True, text=True)
            if result.returncode == 0:
                r_home = result.stdout.strip()
        except:
            pass
    
    if r_home:
        os.environ['R_HOME'] = r_home
        print(f"[INFO] R_HOME set to: {r_home}")
    
    # Set R_USER
    os.environ['R_USER'] = os.path.expanduser('~')
    
    # Set R library path
    r_libs_user = os.path.expanduser('~/R/library')
    os.environ['R_LIBS_USER'] = r_libs_user
    os.makedirs(r_libs_user, exist_ok=True)


def mclust_R_fixed(adata, num_cluster, used_obsm='X_pca'):
    """
    Fixed version of mclust_R function for clustering.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object
    num_cluster : int
        Number of clusters
    used_obsm : str
        Key in obsm to use for clustering
        
    Returns
    -------
    anndata.AnnData
        Updated adata with clustering results
    """
    print(f"[INFO] Running mclust clustering with {num_cluster} clusters...")
    
    # Get the data
    if used_obsm in adata.obsm:
        cluster_data = adata.obsm[used_obsm]
    else:
        print(f"[WARNING] {used_obsm} not found, using SpatialGlue embedding")
        cluster_data = adata.obsm['SpatialGlue']
    
    # Convert to R format using modern rpy2 syntax
    with localconverter(robjects.default_converter + numpy2ri.converter + pandas2ri.converter):
        r_data = robjects.conversion.py2rpy(cluster_data)
    
    # Load mclust in R
    try:
        robjects.r('library(mclust)')
    except Exception as e:
        print(f"[ERROR] Failed to load mclust R package: {e}")
        print("[INFO] Falling back to Leiden clustering...")
        sc.pp.neighbors(adata, use_rep='SpatialGlue', n_neighbors=30)
        sc.tl.leiden(adata, resolution=0.5, key_added='mclust')
        return adata
    
    # Run mclust clustering
    robjects.globalenv['cluster_data'] = r_data
    robjects.globalenv['n_clusters'] = num_cluster
    
    r_script = f"""
    set.seed(2020)
    result <- Mclust(cluster_data, G={num_cluster})
    cluster_labels <- result$classification
    """
    
    try:
        robjects.r(r_script)
        
        # Get results back
        with localconverter(robjects.default_converter + numpy2ri.converter):
            cluster_labels = robjects.conversion.rpy2py(robjects.r['cluster_labels'])
        
        # Add to adata
        adata.obs['mclust'] = cluster_labels.astype(str)
        print(f"[SUCCESS] mclust clustering completed with {len(np.unique(cluster_labels))} clusters")
        
    except Exception as e:
        print(f"[ERROR] mclust clustering failed: {e}")
        print("[INFO] Falling back to Leiden clustering...")
        sc.pp.neighbors(adata, use_rep='SpatialGlue', n_neighbors=30)
        sc.tl.leiden(adata, resolution=0.5, key_added='mclust')
    
    return adata


def preprocess_for_spatialglue(adata_rna, adata_adt):
    """
    Preprocess data for SpatialGlue integration.
    
    Parameters
    ----------
    adata_rna : anndata.AnnData
        RNA data
    adata_adt : anndata.AnnData
        ADT/protein data
        
    Returns
    -------
    tuple
        Preprocessed (adata_rna, adata_adt)
    """
    print("[INFO] Preprocessing data for SpatialGlue...")
    
    # Preprocess RNA data
    print("[INFO] Preprocessing RNA data...")
    sc.pp.filter_genes(adata_rna, min_cells=10)
    sc.pp.highly_variable_genes(adata_rna, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata_rna, target_sum=1e4)
    sc.pp.log1p(adata_rna)
    sc.pp.scale(adata_rna)
    
    adata_rna_high = adata_rna[:, adata_rna.var['highly_variable']]
    adata_rna.obsm['feat'] = pca(adata_rna_high, n_comps=50)
    
    # Preprocess protein data
    print("[INFO] Preprocessing protein data...")
    if scipy.sparse.issparse(adata_adt.X):
        adata_adt.X = adata_adt.X.toarray()
    
    adata_adt = clr_normalize_each_cell(adata_adt)
    sc.pp.scale(adata_adt)
    adata_adt.obsm['feat'] = pca(adata_adt, n_comps=50)
    
    return adata_rna, adata_adt


def run_spatialglue_integration(adata_rna, adata_adt, random_seed=2025):
    """
    Run SpatialGlue integration.
    
    Parameters
    ----------
    adata_rna : anndata.AnnData
        Preprocessed RNA data
    adata_adt : anndata.AnnData
        Preprocessed ADT data
    random_seed : int
        Random seed for reproducibility
        
    Returns
    -------
    anndata.AnnData
        Integrated data with SpatialGlue embeddings
    """
    print("[INFO] Running SpatialGlue integration...")
    
    # Set random seed
    fix_seed(random_seed)
    
    # Set device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[INFO] Using device: {device}")
    
    # Specify data type
    data_type = 'Macsima_Xenium'
    
    # Construct neighbor graph
    print("[INFO] Constructing neighbor graph...")
    data = construct_neighbor_graph(adata_rna, adata_adt, datatype=data_type)
    
    # Train model
    print("[INFO] Training SpatialGlue model...")
    model = Train_SpatialGlue(data, datatype=data_type, device=device)
    output = model.train()
    
    # Create integrated AnnData object
    adata_integrated = adata_rna.copy()
    adata_integrated.obsm['emb_latent_omics1'] = output['emb_latent_omics1']
    adata_integrated.obsm['emb_latent_omics2'] = output['emb_latent_omics2']
    adata_integrated.obsm['SpatialGlue_Mics_Xenium'] = output['SpatialGlue']
    adata_integrated.obsm['SpatialGlue'] = output['SpatialGlue']
    adata_integrated.obsm['alpha'] = output['alpha']
    adata_integrated.obsm['alpha_omics1'] = output['alpha_omics1']
    adata_integrated.obsm['alpha_omics2'] = output['alpha_omics2']
    
    print("[SUCCESS] SpatialGlue integration completed!")
    return adata_integrated


def perform_clustering(adata, n_clusters=6, clustering_method='mclust'):
    """
    Perform clustering on integrated data.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Integrated data
    n_clusters : int
        Number of clusters
    clustering_method : str
        Clustering method ('mclust', 'leiden', 'louvain')
        
    Returns
    -------
    anndata.AnnData
        Data with clustering results
    """
    print(f"[INFO] Performing {clustering_method} clustering with {n_clusters} clusters...")
    
    if clustering_method == 'mclust':
        # Monkey patch the function temporarily
        import SpatialGlue.utils as sg_utils
        sg_utils.mclust_R = mclust_R_fixed
        
        # Run clustering
        clustering(adata, key='SpatialGlue_Mics_Xenium', add_key='SpatialGlue_Mics_Xenium', 
                  method=clustering_method, n_clusters=n_clusters, use_pca=True)
    else:
        # Use scanpy clustering methods
        sc.pp.neighbors(adata, use_rep='SpatialGlue_Mics_Xenium', n_neighbors=30)
        if clustering_method == 'leiden':
            sc.tl.leiden(adata, resolution=0.5, key_added='SpatialGlue_Mics_Xenium')
        elif clustering_method == 'louvain':
            sc.tl.louvain(adata, resolution=0.5, key_added='SpatialGlue_Mics_Xenium')
    
    return adata


def create_visualizations(adata, n_clusters, output_dir):
    """
    Create visualization plots for SpatialGlue results.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Integrated data with clustering
    n_clusters : int
        Number of clusters
    output_dir : Path
        Output directory for plots
    """
    print("[INFO] Creating visualization plots...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Fix spatial coordinates orientation
    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]
    
    # Compute UMAP
    sc.pp.neighbors(adata, use_rep='SpatialGlue_Mics_Xenium', n_neighbors=30)
    sc.tl.umap(adata)
    
    # 1. Large overview plot
    fig, ax_list = plt.subplots(1, 2, figsize=(30, 15))
    sc.pl.umap(adata, color='SpatialGlue_Mics_Xenium', ax=ax_list[0], 
               title='SpatialGlue_Mics_Xenium', s=20, show=False)
    sc.pl.embedding(adata, basis='spatial', color='SpatialGlue_Mics_Xenium', 
                   ax=ax_list[1], title='SpatialGlue_Mics_Xenium', s=20, show=False)
    
    plt.tight_layout(w_pad=0.3)
    plt.savefig(output_dir / f"V3_MACsimas_vs_Xenium_cluster_{n_clusters}.png", 
                dpi=300, bbox_inches="tight")
    plt.close()
    
    # 2. UMAP only plot
    fig, ax = plt.subplots(figsize=(30, 30))
    sc.pl.umap(adata, color='SpatialGlue_Mics_Xenium', ax=ax, 
               title='UMAP of SpatialGlue_Mics_Xenium', s=30, show=False)
    plt.savefig(output_dir / f"umap_plot_V1_MACsimas_vs_Xenium_n_clusters{n_clusters}.png", 
                dpi=300, bbox_inches="tight")
    plt.close()
    
    # 3. Compact annotation plot
    fig, ax_list = plt.subplots(1, 2, figsize=(7, 5))
    sc.pl.umap(adata, color='SpatialGlue_Mics_Xenium', ax=ax_list[0], 
               title='SpatialGlue_Mics_Xenium', s=15, show=False)
    sc.pl.embedding(adata, basis='spatial', color='SpatialGlue_Mics_Xenium', 
                   ax=ax_list[1], title='SpatialGlue_Mics_Xenium', s=15, show=False)
    
    ax_list[0].get_legend().remove()
    plt.tight_layout(w_pad=0.3)
    plt.savefig(output_dir / f"With_Annotation_SpatialGlue_Mics_Xenium_n_clusters{n_clusters}.png", 
                dpi=300, bbox_inches="tight")
    plt.close()
    
    # 4. Weight distribution plot
    plt.rcParams['figure.figsize'] = (15, 5)
    df = pd.DataFrame(columns=['RNA Xenium', 'Protein MACsima', 'label'])
    df['RNA Xenium'], df['Protein MACsima'] = adata.obsm['alpha'][:, 0], adata.obsm['alpha'][:, 1]
    df['label'] = adata.obs['SpatialGlue_Mics_Xenium'].values
    df = df.set_index('label').stack().reset_index()
    df.columns = ['label_SpatialGlue', 'Modality', 'Weight value']
    
    fig, ax = plt.subplots(figsize=(15, 5))
    sns.violinplot(data=df, x='label_SpatialGlue', y='Weight value', hue="Modality",
                  split=True, inner="quart", linewidth=1, ax=ax)
    ax.set_title('MACsima vs Xenium with cluster ' + str(n_clusters))
    ax.set_xlabel('SpatialGlue label')
    ax.legend(bbox_to_anchor=(1.4, 1.01), loc='upper right')
    
    plt.tight_layout(w_pad=0.05)
    plt.savefig(output_dir / f"Weight_RNA_vs_protein_final_cluster_{n_clusters}.png", 
                dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"[SUCCESS] Visualizations saved to {output_dir}")


def main_spatialglue_pipeline(rna_path, adt_path, output_dir=None, n_clusters=6, 
                             clustering_method='mclust', random_seed=2025):
    """
    Complete SpatialGlue integration pipeline.
    
    Parameters
    ----------
    rna_path : str or Path
        Path to RNA AnnData file
    adt_path : str or Path
        Path to ADT AnnData file
    output_dir : str or Path, optional
        Output directory
    n_clusters : int
        Number of clusters
    clustering_method : str
        Clustering method
    random_seed : int
        Random seed
        
    Returns
    -------
    anndata.AnnData
        Integrated data with clustering results
    """
    print("=" * 60)
    print("SPATIALGLUE INTEGRATION PIPELINE")
    print("=" * 60)
    
    if output_dir is None:
        output_dir = OUTPUT_BASE_DIR
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Setup environment
        setup_environment()
        
        # Load data
        print(f"[INFO] Loading RNA data from: {rna_path}")
        adata_rna = sc.read_h5ad(rna_path)
        
        print(f"[INFO] Loading ADT data from: {adt_path}")
        adata_adt = sc.read_h5ad(adt_path)
        
        # Make variable names unique
        adata_rna.var_names_make_unique()
        adata_adt.var_names_make_unique()
        
        print(f"[INFO] RNA data shape: {adata_rna.shape}")
        print(f"[INFO] ADT data shape: {adata_adt.shape}")
        
        # Preprocess data
        adata_rna, adata_adt = preprocess_for_spatialglue(adata_rna, adata_adt)
        
        # Run SpatialGlue integration
        adata_integrated = run_spatialglue_integration(adata_rna, adata_adt, random_seed)
        
        # Save intermediate results
        intermediate_path = output_dir / "data_Macsim_Xenium_spatialglue_before_clustering.h5ad"
        print(f"[INFO] Saving intermediate results to: {intermediate_path}")
        adata_integrated.write(intermediate_path, compression="gzip")
        
        # Perform clustering
        adata_integrated = perform_clustering(adata_integrated, n_clusters, clustering_method)
        
        # Create visualizations
        create_visualizations(adata_integrated, n_clusters, output_dir / "plots")
        
        # Save final results
        final_path = output_dir / f"data_Macsim_Xenium_spatialglue_cluster_{n_clusters}.h5ad"
        print(f"[INFO] Saving final results to: {final_path}")
        adata_integrated.write(final_path, compression="gzip")
        
        print("=" * 60)
        print("SPATIALGLUE INTEGRATION COMPLETED SUCCESSFULLY!")
        print(f"Final integrated data saved to: {final_path}")
        print(f"Visualizations saved to: {output_dir / 'plots'}")
        print("=" * 60)
        
        return adata_integrated
        
    except Exception as e:
        print(f"[ERROR] SpatialGlue pipeline failed: {str(e)}")
        raise


if __name__ == "__main__":
    # Example usage
    rna_path = OUTPUT_BASE_DIR / "adata_RNA_xenium.h5ad"
    adt_path = OUTPUT_BASE_DIR / "adata_ADT_mics.h5ad"
    
    if rna_path.exists() and adt_path.exists():
        adata_integrated = main_spatialglue_pipeline(rna_path, adt_path)
    else:
        print("[ERROR] Input files not found. Please run preprocessing.py first.")
