import os
import sys
import argparse
from Source_error import import_error as iE
def parse_args():
    parser = argparse.ArgumentParser(description="-o output_directory path,-n number of sequences")
    parser.add_argument(
        "-o","--output",
        type=str,
        default="./output",  # default output directory
        help="the output directory path, default: ./output"
    )
    parser.add_argument(
        "-n","--number",
        type=int,
        help="the number of sequences is needed"
    )
    # ---- SHM parameters (applied to IGH only, read from env if in subprocess) ----
    parser.add_argument("--shm-rate", type=float,
                        default=float(os.environ.get("TRBSEQMOCK_SHM_RATE", 0.05)),
                        help="Global SHM rate per base (default: 0.05)")
    parser.add_argument("--shm-hotspot-rr", type=float,
                        default=float(os.environ.get("TRBSEQMOCK_SHM_HOTSPOT_RR", 4.0)),
                        help="Hotspot relative risk (default: 4.0)")
    parser.add_argument("--shm-cdr-bias", type=float,
                        default=float(os.environ.get("TRBSEQMOCK_SHM_CDR_BIAS", 1.5)),
                        help="CDR bias vs FR (default: 1.5)")
    parser.add_argument("--shm-cold-spot-rate", type=float,
                        default=float(os.environ.get("TRBSEQMOCK_SHM_COLD_SPOT_RATE", 0.002)),
                        help="Cold spot residual rate (default: 0.002)")
    parser.add_argument("--shm-seed", type=int,
                        default=None,
                        help="SHM random seed")
    return parser.parse_args()

def ensure_output_dir(output_path: str) -> str:
    abs_output_path = os.path.abspath(output_path)
    try:
        os.makedirs(abs_output_path, exist_ok=True)
        print(f"The output will be placed in: {abs_output_path}")
        return abs_output_path
    except PermissionError:
        print(f"Error: No permission to write output to {abs_output_path}", file=sys.stderr)
        sys.exit(1)
    except OSError as e:
        print(f"Failed to create directory {abs_output_path}, reason: {e}", file=sys.stderr)
        sys.exit(1)

def write_file(output_dir: str, filename: str, content: str,) -> None:
    file_path = os.path.join(output_dir, filename)
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"File {filename} has been written to: {file_path}")
    except Exception as e:
        print(f"Error: Failed to write file {file_path}, reason: {e}", file=sys.stderr)

def generate_fasta(v_name, d_name, j_name, sequence, fasta_path_and_name,
                   shm_report: dict = None):
    """
    Generate or append to a standard FASTA file with optional SHM annotation.
    """
    if shm_report is None:
        shm_report = {}

    seq_length = len(sequence)
    shm_count = shm_report.get("count", 0)
    v_identity = shm_report.get("v_identity", 100.0)
    shm_rate = shm_report.get("shm_rate_used", 0.0)

    # Enhanced header: >V|D|J|len|shm_count|v_identity|shm_rate
    header = (f">{v_name}|{d_name}|{j_name}|{seq_length}"
              f"|{shm_count}|{v_identity}|{shm_rate}")
    formatted_seq = sequence

    fasta_entry = f"\n{header}\n{formatted_seq}"

    with open(fasta_path_and_name, 'a+', encoding='utf-8') as f:
        f.seek(0)
        if not f.read(1):
            f.write(header + '\n' + formatted_seq)
        else:
            f.write(fasta_entry)

def write_detail_tsv(
        tsv_path: str,
        v_name: str,
        d_name: str,
        j_name: str,
        v3_del_len: int,
        d5_del_len: int,
        d3_del_len: int,
        j5_del_len: int,
        vd_insert_len: int,
        dj_insert_len: int,
        vd_insert_seq: str,
        dj_insert_seq: str,
        shm_report: dict = None,
        germline_v_seq: str = "",
        germline_j_seq: str = "",
):
    if shm_report is None:
        shm_report = {}

    is_first_line = not os.path.exists(tsv_path)
    with open(tsv_path, "a", encoding="utf-8") as f:
        if is_first_line:
            header = [
                "V_gene", "D_gene", "J_gene",
                "V3_del", "D5_del", "D3_del", "J5_del",
                "VD_ins_len", "DJ_ins_len", "VD_ins_seq", "DJ_ins_seq",
                "SHM_count", "SHM_v_identity", "SHM_cdr_count", "SHM_fr_count",
                "SHM_hotspot_count", "SHM_transitions", "SHM_transversions",
                "SHM_rate_used", "SHM_positions", "SHM_details",
                "Germline_V_seq", "Germline_J_seq",
            ]
            f.write("\t".join(header) + "\n")

        c = shm_report.get("count", 0)
        row_data = [
            str(v_name), str(d_name), str(j_name),
            str(v3_del_len), str(d5_del_len), str(d3_del_len), str(j5_del_len),
            str(vd_insert_len), str(dj_insert_len), vd_insert_seq, dj_insert_seq,
            str(c),
            str(shm_report.get("v_identity", 100.0)),
            str(shm_report.get("cdr_count", 0)),
            str(shm_report.get("fr_count", 0)),
            str(shm_report.get("hotspot_count", 0)),
            str(shm_report.get("transitions", 0)),
            str(shm_report.get("transversions", 0)),
            str(shm_report.get("shm_rate_used", 0.0)),
            ",".join(str(p) for p in shm_report.get("positions", [])),
            ";".join(shm_report.get("details", [])),
            germline_v_seq,
            germline_j_seq,
        ]
        f.write("\t".join(row_data) + "\n")

insertion_stats = {}
deletion_stats = {}

def update_insertion_stats(insertion_type: str, length: int):
    if insertion_type not in insertion_stats:
        insertion_stats[insertion_type] = {}
    insertion_stats[insertion_type][length] = insertion_stats[insertion_type].get(length, 0) + 1

def update_deletion_stats(deletion_type: str, length: int):
    if deletion_type not in deletion_stats:
        deletion_stats[deletion_type] = {}
    deletion_stats[deletion_type][length] = deletion_stats[deletion_type].get(length, 0) + 1

def write_distribution_file(file_path: str, stats_data: dict):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("Type\tLength\tFrequency\tPercentage\n")
        for type_name, length_counts in stats_data.items():
            type_total = sum(length_counts.values())
            for length in sorted(length_counts.keys()):
                freq = length_counts[length]
                percentage = (freq / type_total) * 100 if type_total != 0 else 0.0
                f.write(f"{type_name}\t{length}\t{freq}\t{percentage:.6f}\n")


# =========================================================================
# SHM Statistics accumulators
# =========================================================================
shm_stats_accumulator = {}   # per V-gene: {count_sum, identity_sum, cdr_sum, fr_sum, ...}
shm_per_position = {}        # per V-gene: {IMGT_pos: mutation_count}


def update_shm_stats(v_gene: str, shm_report: dict):
    """Accumulate SHM statistics per V-gene for the summary report."""
    if v_gene not in shm_stats_accumulator:
        shm_stats_accumulator[v_gene] = {
            "seq_count": 0,
            "mut_count_sum": 0,
            "cdr_sum": 0,
            "fr_sum": 0,
            "hotspot_sum": 0,
            "identity_sum": 0.0,
        }
    s = shm_stats_accumulator[v_gene]
    s["seq_count"] += 1
    s["mut_count_sum"] += shm_report.get("count", 0)
    s["cdr_sum"] += shm_report.get("cdr_count", 0)
    s["fr_sum"] += shm_report.get("fr_count", 0)
    s["hotspot_sum"] += shm_report.get("hotspot_count", 0)
    s["identity_sum"] += shm_report.get("v_identity", 100.0)

    # Per-position accumulation
    if v_gene not in shm_per_position:
        shm_per_position[v_gene] = {}
    for pos in shm_report.get("positions", []):
        shm_per_position[v_gene][pos] = shm_per_position[v_gene].get(pos, 0) + 1


def write_shm_statistics_file(file_path: str, stats: dict):
    """Write per-V-gene SHM statistics."""
    if not stats:
        return
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("\t".join([
            "V_gene", "Seq_count",
            "Avg_mutations", "Avg_v_identity",
            "Avg_CDR_mut", "Avg_FR_mut",
            "CDR_FR_ratio", "Hotspot_enrichment",
        ]) + "\n")
        for gene, s in sorted(stats.items()):
            n = s["seq_count"]
            avg_mut = s["mut_count_sum"] / n if n else 0
            avg_id = s["identity_sum"] / n if n else 100
            avg_cdr = s["cdr_sum"] / n if n else 0
            avg_fr = s["fr_sum"] / n if n else 0
            cdr_fr_ratio = avg_cdr / avg_fr if avg_fr > 0 else float('nan')
            hotspot_enrich = s["hotspot_sum"] / s["mut_count_sum"] if s["mut_count_sum"] > 0 else 0
            f.write("\t".join(map(str, [
                gene, n,
                round(avg_mut, 2), round(avg_id, 2),
                round(avg_cdr, 2), round(avg_fr, 2),
                round(cdr_fr_ratio, 3), round(hotspot_enrich, 3),
            ])) + "\n")


# =========================================================================
# Ground Truth JSON writer
# =========================================================================
import json

def write_ground_truth(
    json_path: str,
    seq_id: int,
    v_name: str,
    d_name: str,
    j_name: str,
    germline_v_seq: str,
    germline_j_seq: str,
    v3_del_len: int,
    d5_del_len: int,
    d3_del_len: int,
    j5_del_len: int,
    vd_insert_len: int,
    dj_insert_len: int,
    vd_insert_seq: str,
    dj_insert_seq: str,
    shm_report: dict,
):
    """Write a ground truth record to a JSON Lines file (one record per line)."""
    record = {
        "seq_id": seq_id,
        "germline_v_gene": v_name,
        "germline_d_gene": d_name,
        "germline_j_gene": j_name,
        "germline_v_seq": germline_v_seq,
        "germline_j_seq": germline_j_seq,
        "shm_applied": shm_report.get("count", 0) > 0,
        "shm_rate_used": shm_report.get("shm_rate_used", 0.0),
        "shm_count": shm_report.get("count", 0),
        "shm_v_identity": shm_report.get("v_identity", 100.0),
        "rearrangement": {
            "v3_del": v3_del_len,
            "d5_del": d5_del_len,
            "d3_del": d3_del_len,
            "j5_del": j5_del_len,
            "vd_insert_len": vd_insert_len,
            "dj_insert_len": dj_insert_len,
            "vd_insert_seq": vd_insert_seq,
            "dj_insert_seq": dj_insert_seq,
        },
    }
    with open(json_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")