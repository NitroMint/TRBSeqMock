import random

def select_lottory(gene_pool):
    lottory_number = random.randint(1, len(gene_pool))-1
    name_seq_type_bonus=gene_pool[lottory_number]
    return name_seq_type_bonus
