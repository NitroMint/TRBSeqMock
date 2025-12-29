import os
import sys
import argparse

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

def generate_fasta(v_name, d_name, j_name, sequence, fasta_path_and_name):
    """
    Generate or append to a standard FASTA file
    - Create file if it does not exist
    - Append new entry if file exists
    - Sequence is displayed in full on one line without automatic line breaks
    """
    seq_length = len(sequence)
    # Construct FASTA header (unchanged)
    header = f">{v_name}|{d_name}|{j_name}|{seq_length}"
    # Key modification: Remove line break logic, use original sequence directly (no splitting)
    formatted_seq = sequence  # No longer split by 60 characters, sequence remains in one line
    # Combine header and sequence (keep blank line separator between entries)
    fasta_entry = f"\n{header}\n{formatted_seq}"

    # Write/append file logic (unchanged)
    with open(fasta_path_and_name, 'a+', encoding='utf-8') as f:
        # Check if file is empty (avoid blank line before first entry)
        f.seek(0)
        if not f.read(1):  # File is empty
            f.write(header + '\n' + formatted_seq)
        else:  # File has content, append new entry (with blank line separator)
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
        dj_insert_seq: str
):
    is_first_line = not os.path.exists(tsv_path)
    with open(tsv_path, "a", encoding="utf-8") as f:
        if is_first_line:
            header = [
                "V_gene", "D_gene", "J_gene",
                "V3_del", "D5_del", "D3_del", "J5_del",
                "VD_ins", "DJ_ins", "VD_ins", "DJ_ins"
            ]
            f.write("\t".join(header) + "\n")

        row_data = [
            str(v_name), str(d_name), str(j_name),
            str(v3_del_len), str(d5_del_len), str(d3_del_len), str(j5_del_len),
            str(vd_insert_len), str(dj_insert_len), vd_insert_seq, dj_insert_seq
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