# Source_SHM - Somatic Hypermutation Simulation Module
#
# Provides a literature-grounded, position-specific SHM model for IGHV genes.
#
# Key references:
#   Rogozin & Diaz, J Immunol 2004;172(6):3382-4  (DGYW/WRCH hotspot motifs)
#   Di Noia & Neuberger, Annu Rev Biochem 2007;76:1-22  (SHM mechanism review)
#   Kleinstein et al., J Immunol 2003;171(9):4639-49  (mutation rate estimates)
#   Sheng et al., Front Immunol 2017;8:537  (gene-specific substitution profiles)
#   Kirik et al., Front Immunol 2017;8:1433  (germline gene mutation paths)
#   Schramm & Douek, Front Immunol 2018;9:1876  (cold spots, non-hotspot mutations)

from .shm_simulator import apply_shm, SHMConfig, scan_hotspots, compute_hotspot_enrichment
from .shm_constants import (
    WRC_PATTERN, GYW_PATTERN,
    IMGT_CDR_FR_BOUNDARIES,
    get_region_for_position,
    is_cdr_position,
    DEFAULT_SHM_RATE,
    DEFAULT_HOTSPOT_RR,
    DEFAULT_CDR_BIAS,
    DEFAULT_COLD_SPOT_RATE,
)
