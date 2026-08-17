"""Python API for the C++ cellAdmix core."""

from .dataset import CellAdmix
from .factor_sources import CellAdmixFactorSourceScore
from .fit import CellAdmixFit
from .score import CellAdmixScore
from .correction import CellAdmixCorrection
from .audit import CellAdmixAudit, CellAdmixCleanupReport
from .io import read_annotation
from .plotting import plot_spatial
from .spatialdata import add_corrected_counts_to_spatialdata, add_fit_to_spatialdata, from_spatialdata
from .state import annotation_knn_purity

__version__ = "0.0.1"

__all__ = [
    "CellAdmix",
    "CellAdmixFactorSourceScore",
    "CellAdmixFit",
    "CellAdmixScore",
    "CellAdmixCorrection",
    "read_annotation",
    "plot_spatial",
    "annotation_knn_purity",
    "from_spatialdata",
    "add_fit_to_spatialdata",
    "add_corrected_counts_to_spatialdata",
]
