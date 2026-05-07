from .read_fasta import v_pool,d_pool,j_pool
from .read_distributions import v3_dist,dj_dist

"""
note: the co-structure of v_pool,d_pool,j_pool be like  (gene_names[idx], sequences[idx], genotypes[idx])
It differs from the version in TRBSeqMock which only has three D gene and at that situation I ignored the genotype difference.

v3_dist ,dj_dist be like ()
"""