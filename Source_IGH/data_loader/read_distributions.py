from .config import DELETION_DIST_FILE, INSERTION_DIST_FILE

"""
In this script, we got the P-L distribution of the deletion and insertion
in our files in the direction Ref
We choose the V3's deletion and DJ's insertion
"""

def read_distribution_files():
    # Read deletion distribution file and match V3 data by gene name (instead of hardcoding index)
    with open(DELETION_DIST_FILE, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
        v3_dist = None
        for line in lines:
            if line.startswith('V3'):  # Match by gene name to avoid index error
                v3_row = line.split()
                v3_dist = list(map(float, v3_row[1:]))  # Skip gene name, extract numeric values
                break
        if v3_dist is None:
            raise ValueError("V3 corresponding row not found in deletion distribution file")


    # Read insertion distribution file and match DJ data by gene name

    with open(INSERTION_DIST_FILE, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
        dj_dist = None
        for line in lines:
            if line.startswith('DJ'):  # Match by gene name for better robustness
                dj_row = line.split()
                dj_dist = list(map(float, dj_row[1:]))  # Skip gene name, extract numeric values
                break
        if dj_dist is None:
            raise ValueError("DJ corresponding row not found in insertion distribution file")

    return v3_dist, dj_dist

v3_dist, dj_dist = read_distribution_files()
