from seq_sim import main, generate_fasta, ensure_output_dir, parse_args,write_detail_tsv,update_deletion_stats,update_insertion_stats,write_distribution_file,insertion_stats,deletion_stats
import os

args = parse_args()
output_dir = args.output
num_repeats = args.number

ensure_output_dir(output_dir)

fasta_path = os.path.join(output_dir, "result.fasta")
tsv_path = os.path.join(output_dir, "detail_result.tsv")
insert_dist_path = os.path.join(output_dir, "insertion_distribution.tsv")
delete_dist_path = os.path.join(output_dir, "deletion_distribution.tsv")

for repeat in range(num_repeats):
    (
        v, d, j, seq,
        v3_del_len, d5_del_len, d3_del_len, j5_del_len,
        vd_insert_len, dj_insert_len, vd_insert_seq, dj_insert_seq
    ) = main()


    generate_fasta(v, d, j, seq, fasta_path)


    write_detail_tsv(
        tsv_path=tsv_path,
        v_name=v, d_name=d, j_name=j,
        v3_del_len=v3_del_len, d5_del_len=d5_del_len, d3_del_len=d3_del_len, j5_del_len=j5_del_len,
        vd_insert_len=vd_insert_len, dj_insert_len=dj_insert_len,
        vd_insert_seq=vd_insert_seq, dj_insert_seq=dj_insert_seq
    )


    update_insertion_stats("VD_Insertion", vd_insert_len)
    update_insertion_stats("DJ_Insertion", dj_insert_len)


    update_deletion_stats("V3_Del", v3_del_len)
    update_deletion_stats("D5_Del", d5_del_len)
    update_deletion_stats("D3_Del", d3_del_len)
    update_deletion_stats("J5_Del", j5_del_len)

write_distribution_file(insert_dist_path, insertion_stats)
write_distribution_file(delete_dist_path, deletion_stats)

print(f"Success! Generated {num_repeats} sequences to {os.path.abspath(fasta_path)}")
print(f"Detail info saved to {os.path.abspath(tsv_path)}")
print(f"Insertion distribution saved to {os.path.abspath(insert_dist_path)}")
print(f"Deletion distribution saved to {os.path.abspath(delete_dist_path)}")