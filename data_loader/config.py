import os

# Get package directory for path alignment
PROCESS_DIR = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(os.path.dirname(PROCESS_DIR), "Ref")


# Distribution file paths
DELETION_DIST_FILE = os.path.join(REF, "VDJ_deletion_nt_len.dis.txt")
INSERTION_DIST_FILE = os.path.join(REF, "VDJ_insertion_nt_len.dis.txt")

#VDJ fasta file paths
V_FASTA = os.path.join(REF, "TRBV.fasta.txt")
D_FASTA = os.path.join(REF, "TRBD.fasta.txt")
J_FASTA = os.path.join(REF, "TRBJ.fasta.V2.txt")
