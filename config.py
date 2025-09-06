"""
Configuration file for MACsima and Xenium data integration pipeline.
Modify paths and parameters according to your data setup.
"""

from pathlib import Path

# ============================================================================
# DATA PATHS
# ============================================================================
MACSIMA_DIR = Path("/home/fenosoa/links/scratch/Data_integration/MACSimaROI1/ROI1/_clean_min")
XENIUM_DIR = Path("/home/fenosoa/links/scratch/Data_integration/XeniumROI1/output-XETG00325__0029986__Region_1__20241108__233229")
OUTPUT_BASE_DIR = Path("/home/fenosoa/links/scratch/Data_integration/MACSimaROI1/ROI1")

# ============================================================================
# OUTPUT PATHS
# ============================================================================
OUT_ZARR = OUTPUT_BASE_DIR / "macsima_aggregated.zarr"
XENIUM_EXPLORER_OUTPUT = OUTPUT_BASE_DIR / "xenium_explorer_output.explorer"
XENIUM_EXPLORER_SPATIALGLUE = OUTPUT_BASE_DIR / "xenium_explorer_output_MACSima_Xenium_SpatialGlue.explorer"

# ============================================================================
# SOPA CONFIGURATION
# ============================================================================
IMAGE_KEY = "_clean_min"
CELLS_KEY = "CELLPOSE_BOUNDARIES"
GENES_TABLE_KEY = "table_genes"
PROTEINS_TABLE_KEY = "table_proteins"
COMBINED_TABLE_KEY = "combined"

# ============================================================================
# CELLPOSE SEGMENTATION PARAMETERS
# ============================================================================
DIAMETER_PX = 18
RUN_WITH_GPU = True
PRETRAINED_MODEL = "/home/fenosoa/.cellpose/models/cpsam"
PATCH_WIDTH = 2000
PATCH_OVERLAP = 50

# ============================================================================
# PROCESSING SETTINGS
# ============================================================================
import sopa
sopa.settings.parallelization_backend = 'dask'
