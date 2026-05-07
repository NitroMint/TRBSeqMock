from seq_sim import main, SHM_PARAMS
import os
from Source_write import*



args = parse_args()
output_dir = args.output
num_repeats = args.number

# Inject SHM parameters from CLI args into the main module
SHM_PARAMS["shm_rate"] = getattr(args, "shm_rate", 0.05)
SHM_PARAMS["hotspot_rr"] = getattr(args, "shm_hotspot_rr", 4.0)
SHM_PARAMS["cdr_bias"] = getattr(args, "shm_cdr_bias", 1.5)
SHM_PARAMS["cold_spot_rate"] = getattr(args, "shm_cold_spot_rate", 0.002)
SHM_PARAMS["seed"] = getattr(args, "shm_seed", None)

ensure_output_dir(output_dir)

fasta_path = os.path.join(output_dir, "result.fasta")
tsv_path = os.path.join(output_dir, "detail_result.tsv")
insert_dist_path = os.path.join(output_dir, "insertion_distribution.tsv")
delete_dist_path = os.path.join(output_dir, "deletion_distribution.tsv")
shm_stat_path = os.path.join(output_dir, "shm_statistics.tsv")
ground_truth_path = os.path.join(output_dir, "ground_truth.json")

for repeat in range(num_repeats):
    (
        v, d, j, seq,
        v3_del_len, d5_del_len, d3_del_len, j5_del_len,
        vd_insert_len, dj_insert_len, vd_insert_seq, dj_insert_seq,
        shm_report,
        germline_v_seq,
        germline_j_seq,
    ) = main()


    generate_fasta(v, d, j, seq, fasta_path, shm_report)


    write_detail_tsv(
        tsv_path=tsv_path,
        v_name=v, d_name=d, j_name=j,
        v3_del_len=v3_del_len, d5_del_len=d5_del_len, d3_del_len=d3_del_len, j5_del_len=j5_del_len,
        vd_insert_len=vd_insert_len, dj_insert_len=dj_insert_len,
        vd_insert_seq=vd_insert_seq, dj_insert_seq=dj_insert_seq,
        shm_report=shm_report,
        germline_v_seq=germline_v_seq,
        germline_j_seq=germline_j_seq,
    )


    update_insertion_stats("VD_Insertion", vd_insert_len)
    update_insertion_stats("DJ_Insertion", dj_insert_len)


    update_deletion_stats("V3_Del", v3_del_len)
    update_deletion_stats("D5_Del", d5_del_len)
    update_deletion_stats("D3_Del", d3_del_len)
    update_deletion_stats("J5_Del", j5_del_len)

    update_shm_stats(v, shm_report)

    write_ground_truth(
        ground_truth_path, repeat,
        v, d, j,
        germline_v_seq, germline_j_seq,
        v3_del_len, d5_del_len, d3_del_len, j5_del_len,
        vd_insert_len, dj_insert_len, vd_insert_seq, dj_insert_seq,
        shm_report,
    )


write_distribution_file(insert_dist_path, insertion_stats)
write_distribution_file(delete_dist_path, deletion_stats)
write_shm_statistics_file(shm_stat_path, shm_stats_accumulator)

print(f"Success! Generated {num_repeats} sequences to {os.path.abspath(fasta_path)}")
print(f"Detail info saved to {os.path.abspath(tsv_path)}")
print(f"SHM statistics saved to {os.path.abspath(shm_stat_path)}")
print(f"Ground truth saved to {os.path.abspath(ground_truth_path)}")
print(f"Insertion distribution saved to {os.path.abspath(insert_dist_path)}")
print(f"Deletion distribution saved to {os.path.abspath(delete_dist_path)}")