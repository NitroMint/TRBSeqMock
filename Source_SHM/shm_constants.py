"""
SHM biological constants and motif definitions.

All constants are grounded in published literature (see module docstrings for
specific references).  IMGT unique numbering for IGHV is used for CDR/FR
boundary definition.
"""

import re

# ---------------------------------------------------------------------------
# Hotspot motifs (Rogozin & Diaz, J Immunol 2004)
#   WRC  (W = A/T, R = A/G)  -- AID targets the C in this context
#   GYW  (Y = C/T, W = A/T)  -- reverse complement; G is the hotspot
# The canonical extended form DGYW/WRCH is more predictive for G:C mutations,
# but WRC/GYW captures the core motif with adequate sensitivity for simulation.
# ---------------------------------------------------------------------------
WRC_PATTERN = re.compile(r"[AT][AG]C")   # W = [AT], R = [AG], core=C
GYW_PATTERN = re.compile(r"G[CT][AT]")   # core=G, Y = [CT], W = [AT]

# Extended DGYW/WRCH patterns (used for finer-grain risk)
DGYW_PATTERN = re.compile(r"[AGT]G[CT][AT]")    # D = not C
WRCH_PATTERN = re.compile(r"[AT][AG]C[ACT]")    # H = not G


# ---------------------------------------------------------------------------
# IMGT unique numbering boundaries for IGHV CDR and FR regions.
# Positions use the IMGT unique numbering system (Lefranc et al., 2003).
# FR1: 1-26, CDR1: 27-38, FR2: 39-55, CDR2: 56-65, FR3: 66-104
# Beyond position 104 is CDR3 (excluded from SHM target).
# ---------------------------------------------------------------------------
IMGT_CDR_FR_BOUNDARIES = {
    "FR1":  (1, 26),
    "CDR1": (27, 38),
    "FR2":  (39, 55),
    "CDR2": (56, 65),
    "FR3":  (66, 104),
}

# IMGT positions that are considered CDR (hotter regions for mutation)
CDR_POSITIONS = set()
for start, end in [(27, 38), (56, 65)]:
    for p in range(start, end + 1):
        CDR_POSITIONS.add(p)

FR_POSITIONS = set()
for start, end in [(1, 26), (39, 55), (66, 104)]:
    for p in range(start, end + 1):
        FR_POSITIONS.add(p)

# All IGHV target positions (FR1-CDR1-FR2-CDR2-FR3)
IGHV_TARGET_POSITIONS = CDR_POSITIONS | FR_POSITIONS


# ---------------------------------------------------------------------------
# Base substitution probability weights for SHM.
#
# AID (Activation-Induced Deaminase) deaminates C → U on single-stranded DNA.
# Repair of the U:G mismatch can yield:
#   - C→T (transition, uracil replaced by T via replication)
#   - C→G or C→A (transversions, via UNG/APE-mediated repair)
# The complementary strand yields G→A, G→C, G→T.
#
# Empirically, transitions account for ~55-60% of SHM events
# (Di Noia & Neuberger, 2007; Sheng et al., 2017).
#
# Weights below are relative probabilities.
# Each row sums to 1.0, normalized later if needed.
# ---------------------------------------------------------------------------

# Substitution weight matrix: [from_base][to_base_index]
# Order: A=0, C=1, G=2, T=3
# Derived from Sheng et al. (2017) aggregate substitution profiles.
SUBSTITUTION_WEIGHTS = {
    "A": {"C": 0.30, "G": 0.35, "T": 0.35},   # A→G transition slightly favored
    "C": {"A": 0.15, "G": 0.10, "T": 0.75},   # C→T transition strongly favored (AID deamination)
    "G": {"A": 0.55, "C": 0.20, "T": 0.25},   # G→A transition favored (complementary C→T)
    "T": {"A": 0.30, "C": 0.35, "G": 0.35},   # roughly balanced
}

# Pre-computed transition indices for quick lookup
TRANSITIONS = {
    "A": "G", "G": "A",
    "C": "T", "T": "C",
}

ALL_BASES = ["A", "C", "G", "T"]


def is_transition(from_base: str, to_base: str) -> bool:
    """Return True if the mutation is a transition (purine↔purine or pyrimidine↔pyrimidine).

    In SHM biochemistry, transitions (C↔T on one strand, G↔A on the complementary
    strand) are favored due to AID-mediated deamination and subsequent repair pathways
    (Di Noia & Neuberger, 2007).
    """
    return to_base == TRANSITIONS.get(from_base)


# ---------------------------------------------------------------------------
# Default SHM configuration values (literature-grounded).
# ---------------------------------------------------------------------------
DEFAULT_SHM_RATE = 0.05             # 5% per base (memory B cell typical, Zheng 2005)
DEFAULT_HOTSPOT_RR = 4.0            # risk ratio for WRC/GYW hotspots vs background
DEFAULT_CDR_BIAS = 1.5              # CDR mutation rate multiplier vs FR
DEFAULT_COLD_SPOT_RATE = 0.002      # residual mutation rate outside hotspots (Schramm 2018)


def get_region_for_position(imgt_pos: int) -> str:
    """Return the region name ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3') for an IMGT position.

    Args:
        imgt_pos: IMGT unique numbering position (1-based)

    Returns:
        Region name string, or 'unknown' if out of IGHV bounds.
    """
    for region, (start, end) in IMGT_CDR_FR_BOUNDARIES.items():
        if start <= imgt_pos <= end:
            return region
    return "unknown"


def is_cdr_position(imgt_pos: int) -> bool:
    """Return True if the IMGT position falls within CDR1 or CDR2.

    Args:
        imgt_pos: IMGT unique numbering position (1-based)

    Returns:
        bool
    """
    return imgt_pos in CDR_POSITIONS
