import random

def normalize_prob(prob_list):
    total = sum(prob_list)
    if total == 0:
        return [1/len(prob_list)] * len(prob_list)
    normalized = [p / total for p in prob_list]
    return normalized

def del_lottory(lp_list_normalized):

    rand_num = random.random()
    cumulative_prob = 0.0
    for idx, prob in enumerate(lp_list_normalized):
        cumulative_prob += prob
        if rand_num <= cumulative_prob:
            return idx
    return 0


def del_process(gene, end, del_len):


    gene_name, seq, genotype = gene
    seq_len = len(seq)
    del_end = end

    if del_len <= 0:
        del_seq = ""
        seq_new = seq
    elif del_len > seq_len:
        del_seq = seq
        seq_new = ''
    else:
        if end == '5':
            del_seq = seq[:del_len]
            seq_new = seq[del_len:]
        else:
            del_seq = seq[-del_len:]
            seq_new = seq[:-del_len]
    real_del_len = len(del_seq)
    gene_new = (gene_name, seq_new, genotype)
    return gene_new, del_end, real_del_len, del_seq
