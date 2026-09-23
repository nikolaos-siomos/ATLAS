"""ATLAS intercomparison configuration and processing helpers."""

from .config import parse_intercomparison_file, parse_intercomparison_ini
from .bundles import prepare_intercomparison_bundles
from .time_processing import filter_and_average_intercomparison_bundles
from .background import apply_intercomparison_background_correction
from .normalization import apply_intercomparison_normalization
from .vertical_processing import harmonize_intercomparison_vertical

__all__ = [
    "parse_intercomparison_file",
    "parse_intercomparison_ini",
    "prepare_intercomparison_bundles",
    "filter_and_average_intercomparison_bundles",
    "apply_intercomparison_background_correction",
    "apply_intercomparison_normalization",
    "harmonize_intercomparison_vertical",
]
