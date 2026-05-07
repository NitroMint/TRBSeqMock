import random

def create_insertion(lp_list_normalized):

    base = ["A", "C", "G", "T"]

    rand_num = random.random()
    cumulative_prob = 0.0
    insert_len = 0
    for idx, prob in enumerate(lp_list_normalized):
        cumulative_prob += prob
        if rand_num <= cumulative_prob:
            insert_len = idx
            break

    if insert_len == 0:
        insert_seq = ""
    else:
        insert_seq = ''.join(random.choice(base) for _ in range(insert_len))

    return insert_seq


