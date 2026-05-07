from .gene_select import select_lottory
import Source_IGH.data_loader as dl
from .delete_process import del_process,del_lottory,normalize_prob
from .create_insertion import create_insertion
from .gene_assemble import assemble
from  Source_error import import_error
from  Source_SHM import apply_shm, SHMConfig

"""
1.choose a group of VDJ
2.process the V gene(3') >v_gene
3.process the D gene(5') >d_gene_mid
4.process the D gene(3') >d_gene
5.process the J gene(5') >j_gene
6.create an VD insertion fragment >vd_seq
7.create an DJ insertion fragment >dj_seq
8 finally assemble the V-VD-D-DJ-J fragments
9.apply somatic hypermutation (SHM) to the IGHV region [IGH only]

SHM configuration is read from the module-level SHM_CONFIG dict, which
can be populated by the caller (e.g. IGHSeqMock.py) before main() is called.
"""

# Module-level SHM configuration — set by IGHSeqMock.py before calling main()
SHM_PARAMS = {
    "shm_rate": 0.05,
    "hotspot_rr": 4.0,
    "cdr_bias": 1.5,
    "cold_spot_rate": 0.002,
    "seed": None,
}


def main():
    #1.select the genes randomly
    v_bonus = select_lottory(dl.v_pool)
    d_bonus = select_lottory(dl.d_pool)
    j_bonus = select_lottory(dl.j_pool)

    #set the const "prob_del"
    prob_del = normalize_prob(dl.v3_dist)
    prob_add = normalize_prob(dl.dj_dist)

    len_del_gene_v3 = del_lottory(prob_del)
    len_del_gene_d5 = del_lottory(prob_del)
    len_del_gene_d3 = del_lottory(prob_del)
    len_del_gene_j5 = del_lottory(prob_del)

    # 2.process the V gene(3') >v_gene
    v_gene,v3_del_end,v3_del_len,v3_del_seq=del_process(v_bonus,"3",len_del_gene_v3)

    # 3.process the D gene(5') >d_gene_mid
    d_gene_mid,d5_del_end,d5_del_len,d5_del_seq=del_process(d_bonus,"5",len_del_gene_d5)

    # 4.process the D gene(3') >d_gene
    d_gene,d3_del_end,d3_del_len,d3_del_seq=del_process(d_gene_mid,"3",len_del_gene_d3)

    # 5.process the J gene(5') >j_gene
    j_gene,j5_del_end,j5_del_len,j5_del_seq=del_process(j_bonus,"5",len_del_gene_j5)

    #6. & 7.create VD and DJ insertion fragment
    vd_seq = create_insertion(prob_add)
    dj_seq = create_insertion(prob_add)

    #8.finally assemble the V-VD-D-DJ-J fragments
    assembled = assemble(v_gene[1], vd_seq, d_gene[1], dj_seq, j_gene[1])

    # ---- 8.5: Apply sequencing error ----
    final_gene_seq = import_error(assembled)

    # ---- 9: Apply SHM to the IGHV region ----
    shm_config = SHMConfig(
        shm_rate=SHM_PARAMS["shm_rate"],
        hotspot_rr=SHM_PARAMS["hotspot_rr"],
        cdr_bias=SHM_PARAMS["cdr_bias"],
        cold_spot_rate=SHM_PARAMS["cold_spot_rate"],
        seed=SHM_PARAMS.get("seed"),
    )
    # Get the germline V sequence (before any deletion)
    germline_v_seq = v_bonus[1]  # v_bonus = (name, seq, genotype)
    final_gene_seq, shm_report = apply_shm(
        final_gene_seq,
        germline_v_seq,
        v3_del_len,
        shm_config,
    )

    return (
        v_gene[0], d_gene[0], j_gene[0], final_gene_seq,
        v3_del_len, d5_del_len, d3_del_len, j5_del_len,
        len(vd_seq), len(dj_seq), vd_seq, dj_seq,
        shm_report,
        germline_v_seq,  # full germline V for ground truth
        j_bonus[1],      # full germline J for ground truth
    )




