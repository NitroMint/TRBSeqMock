import random
from .config import V_FASTA, D_FASTA, J_FASTA
import re

"""
In this script we got the gene pools,and it can process the IMGT standard referencre 
human VDJ fasta,and finally the full gene can be chosen for the following jobs
2025/12/18  
"""


def read_trb_fasta(fasta_path):

    gene_names = []
    sequences = []
    genotypes = []
    full_lines = []

    with open(fasta_path, 'r', encoding='utf-8') as f:
        current_seq = []
        current_gene_name = ""
        current_genotype = ""

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):

                if current_gene_name and current_seq:
                    full_seq = ''.join(current_seq).replace('.', '').replace(' ', '').replace('\n', '').replace('\r', '').replace('-', '').upper()
                    sequences.append(full_seq)
                    gene_names.append(current_gene_name)
                    genotypes.append(current_genotype)
                    current_seq = []

                full_lines.append(line)
                parts = line.split('|')
                current_gene_name = parts[1] if len(parts) >= 2 else ""
                genotype_match = re.search(r'\*(\d+)(\*\d+)?', current_gene_name)
                current_genotype = genotype_match.group(0) if genotype_match else ""
            else:
                current_seq.append(line)

        if current_gene_name and current_seq:
            full_seq = ''.join(current_seq).replace('.', '').replace(' ', '').upper().replace('\n', '').upper()
            sequences.append(full_seq)
            gene_names.append(current_gene_name)
            genotypes.append(current_genotype)

    return gene_names, sequences, genotypes, full_lines

v_names, v_sequences, v_genotypes, v_full_lines = read_trb_fasta(V_FASTA)
d_names, d_sequences, d_genotypes, _ = read_trb_fasta(D_FASTA)
j_names, j_sequences, j_genotypes, j_full_lines = read_trb_fasta(J_FASTA)

def filter_vj(full_lines, gene_names, sequences, genotypes):
    valid_indices = []
    for idx, line in enumerate(full_lines):
        if (("|F|" in line or "|[F]|" in line or "|(F)|" in line) and
            "partial" not in line):
            valid_indices.append(idx)


    valid_sequences = [
        (gene_names[idx], sequences[idx], genotypes[idx])
        for idx in valid_indices
    ]


    pool_01 = []
    pool_other = []
    for gene_info in valid_sequences:
        genotype = gene_info[2]
        if genotype == "*01":
            pool_01.append(gene_info)
        elif genotype.startswith('*02'):
            pool_other.append(gene_info)

    select_num = max(1, round(len(pool_01) * 0.2))
    select_num = min(select_num, len(pool_other))
    pool_02 = random.sample(pool_other, select_num) if pool_other else []

    final_pool = pool_01 + pool_02
    return final_pool

v_pool = filter_vj(v_full_lines, v_names, v_sequences, v_genotypes)
j_pool = filter_vj(j_full_lines, j_names, j_sequences, j_genotypes)
d_pool = [(d_names[i], d_sequences[i].upper(), d_genotypes[i]) for i in range(len(d_names))]