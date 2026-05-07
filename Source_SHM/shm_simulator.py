"""
Somatic Hypermutation (SHM) simulator for IGHV genes.

Literature basis:
  - Rogozin & Diaz (2004) J Immunol 172:3382-4   DGYW/WRCH hotspot motifs
  - Di Noia & Neuberger (2007) Annu Rev Biochem 76:1-22   SHM mechanism
  - Kleinstein et al. (2003) J Immunol 171:4639-49   mutation rate ~1/1000bp/division
  - Sheng et al. (2017) Front Immunol 8:537   gene-specific substitution profiles
  - Schramm & Douek (2018) Front Immunol 9:1876   cold spots, non-canonical targeting

Model: position-specific probability model for IGHV (FR1-CDR1-FR2-CDR2-FR3 only).
SHM is applied to the IGHV region of the assembled sequence; CDR3, D, J and
N-insertions are left unmutated.
"""

import random
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

from .shm_constants import (
    WRC_PATTERN, GYW_PATTERN,
    DGYW_PATTERN, WRCH_PATTERN,
    CDR_POSITIONS, IGHV_TARGET_POSITIONS,
    SUBSTITUTION_WEIGHTS,
    DEFAULT_SHM_RATE,
    DEFAULT_HOTSPOT_RR,
    DEFAULT_CDR_BIAS,
    DEFAULT_COLD_SPOT_RATE,
    get_region_for_position,
    is_transition,
)


# ---------------------------------------------------------------------------
@dataclass
class SHMConfig:
    """SHM simulation parameters.

    Attributes:
        shm_rate: Global per-base mutation probability (default 0.05 = 5%).
        hotspot_rr: Relative risk multiplier for WRC/GYW hotspot positions.
        cdr_bias: Mutation rate multiplier for CDR positions vs FR.
        cold_spot_rate: Residual mutation rate outside hotspots (default 0.002).
        seed: Random seed for reproducibility.  None means system default.
    """
    shm_rate: float = DEFAULT_SHM_RATE
    hotspot_rr: float = DEFAULT_HOTSPOT_RR
    cdr_bias: float = DEFAULT_CDR_BIAS
    cold_spot_rate: float = DEFAULT_COLD_SPOT_RATE
    seed: Optional[int] = None

    def __post_init__(self):
        if self.seed is not None:
            random.seed(self.seed)


# ---------------------------------------------------------------------------
def _find_motif_positions(seq: str, pattern, core_offset: int) -> List[int]:
    """Return 1-based positions of the core hotspot base in a motif."""
    positions = []
    for m in pattern.finditer(seq):
        pos = m.start() + core_offset + 1
        positions.append(pos)
    return positions


def scan_hotspots(v_seq: str) -> Dict[int, str]:
    """Identify WRC/GYW/DGYW/WRCH hotspot positions in a germline IGHV sequence.

    Args:
        v_seq: Uppercase germline IGHV sequence (FR1..FR3, IMGT-aligned).

    Returns:
        Dict mapping 1-based IMGT position -> motif type(s).
    """
    hotspots: Dict[int, str] = {}

    for p in _find_motif_positions(v_seq, WRC_PATTERN, core_offset=2):
        hotspots[p] = "WRC"
    for p in _find_motif_positions(v_seq, GYW_PATTERN, core_offset=0):
        hotspots[p] = hotspots.get(p, "") + "+GYW" if p in hotspots else "GYW"
    for p in _find_motif_positions(v_seq, DGYW_PATTERN, core_offset=1):
        hotspots[p] = hotspots.get(p, "") + "+DGYW" if p in hotspots else "DGYW"
    for p in _find_motif_positions(v_seq, WRCH_PATTERN, core_offset=2):
        hotspots[p] = hotspots.get(p, "") + "+WRCH" if p in hotspots else "WRCH"

    return {p: label for p, label in hotspots.items()
            if p in IGHV_TARGET_POSITIONS}


# ---------------------------------------------------------------------------
def apply_shm(
    assembled_seq: str,
    v_germline_seq: str,
    v3_del_len: int,
    config: Optional[SHMConfig] = None,
) -> Tuple[str, Dict]:
    """Apply SHM to the IGHV portion of an assembled rearranged sequence.

    Only targets FR1-CDR1-FR2-CDR2-FR3 (IMGT positions 1-104).

    Args:
        assembled_seq: Full rearranged sequence (V+VD_ins+D+DJ_ins+J), uppercase.
        v_germline_seq: Germline V gene sequence BEFORE any deletion (uppercase).
        v3_del_len: Bases deleted from V 3' end during recombination.
        config: SHM configuration.  Uses defaults if None.

    Returns:
        (mutated_sequence, shm_report_dict)
    """
    if config is None:
        config = SHMConfig()

    v_region_len = max(0, len(v_germline_seq) - v3_del_len)
    if v_region_len == 0:
        return assembled_seq, _empty_report()

    germline_surviving = v_germline_seq[:v_region_len]
    hotspot_map = scan_hotspots(germline_surviving)

    mutated_seq = list(assembled_seq)
    mutations: List[Dict] = []

    for imgt_pos in range(1, v_region_len + 1):
        seq_idx = imgt_pos - 1
        if seq_idx >= len(assembled_seq):
            break
        original_base = assembled_seq[seq_idx]
        if original_base not in SUBSTITUTION_WEIGHTS:
            continue

        # Per-position probability
        if imgt_pos in hotspot_map:
            p = config.shm_rate * config.hotspot_rr
        else:
            p = config.cold_spot_rate
        if imgt_pos in CDR_POSITIONS:
            p *= config.cdr_bias
        p = min(p, 0.95)

        if random.random() < p:
            weights = SUBSTITUTION_WEIGHTS[original_base]
            bases, w = list(weights.keys()), list(weights.values())
            new_base = random.choices(bases, weights=w, k=1)[0]
            mutated_seq[seq_idx] = new_base
            mutations.append({
                "imgt_pos": imgt_pos,
                "from": original_base, "to": new_base,
                "region": get_region_for_position(imgt_pos),
                "is_hotspot": imgt_pos in hotspot_map,
                "hotspot_motif": hotspot_map.get(imgt_pos, ""),
            })

    cdr_count = sum(1 for m in mutations if m["region"].startswith("CDR"))
    fr_count = sum(1 for m in mutations if m["region"].startswith("FR"))
    hotspot_count = sum(1 for m in mutations if m["is_hotspot"])
    count = len(mutations)
    n_transitions = sum(1 for m in mutations if is_transition(m["from"], m["to"]))
    n_transversions = count - n_transitions

    report = {
        "count": count,
        "cdr_count": cdr_count,
        "fr_count": fr_count,
        "hotspot_count": hotspot_count,
        "transitions": n_transitions,
        "transversions": n_transversions,
        "positions": [m["imgt_pos"] for m in mutations],
        "details": [f"{m['imgt_pos']}:{m['from']}>{m['to']}" for m in mutations],
        "v_identity": round((1 - count / v_region_len) * 100, 2) if v_region_len > 0 else 100.0,
        "v_region_len": v_region_len,
        "shm_rate_used": config.shm_rate,
    }
    return "".join(mutated_seq), report


def _empty_report() -> Dict:
    return {
        "count": 0, "cdr_count": 0, "fr_count": 0, "hotspot_count": 0,
        "transitions": 0, "transversions": 0,
        "positions": [], "details": [],
        "v_identity": 100.0, "v_region_len": 0, "shm_rate_used": 0.0,
    }


def compute_hotspot_enrichment(
    v_germline_seq: str,
    mutation_positions: List[int],
) -> float:
    """Compute hotspot enrichment ratio (observed/expected).

    >1 indicates hotspot enrichment, as expected for real SHM data.
    """
    hotspots = scan_hotspots(v_germline_seq)
    total_positions = min(len(v_germline_seq), max(IGHV_TARGET_POSITIONS, default=104))
    expected_frac = len(hotspots) / total_positions if total_positions > 0 else 0
    if not mutation_positions:
        return 0.0
    observed = sum(1 for p in mutation_positions if p in hotspots) / len(mutation_positions)
    return observed / expected_frac if expected_frac > 0 else float('inf')
