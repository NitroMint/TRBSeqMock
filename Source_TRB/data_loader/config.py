import os.path as path

"""os.path: dirname-father
            join-the dirpath plus another
"""
# Get package directory for path alignment
PROCESS_DIR = path.dirname(path.dirname(path.abspath(__file__)))
Data = path.join(path.dirname(PROCESS_DIR), "Data")

REF = path.join(Data, "Data_Ref")
DIS = path.join(Data, "Data_Distribute")

# Distribution file paths
DELETION_DIST_FILE = path.join(DIS, "TRB_VDJ_deletion_nt_len.dis.txt")
INSERTION_DIST_FILE = path.join(DIS, "TRB_VDJ_insertion_nt_len.dis.txt")

#VDJ fasta file paths
V_FASTA = path.join(REF, "TRBV.fasta.txt")
D_FASTA = path.join(REF, "TRBD.fasta.txt")
J_FASTA = path.join(REF, "TRBJ.fasta.V2.txt")
