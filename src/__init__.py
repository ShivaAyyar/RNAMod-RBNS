"""
Epitranscriptomic RBP Analysis Package

Modules:
    peak_analysis: Peak processing and RBNS scoring
    enrichment_analysis: Modification enrichment testing
    visualization: Plotting functions
    main: Pipeline orchestration
"""

from . import peak_analysis
from . import enrichment_analysis
from . import visualization

__version__ = '1.0.0'
__author__ = 'RNAMod-RBNS Analysis Pipeline'
