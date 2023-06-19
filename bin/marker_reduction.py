import os
import random

from get_LA_info import get_LA_info

def marker_reduction(kmer:int, max_markers:int)->None:
    if not os.path.exists("AFLAP_tmp/LA.txt"):
        raise FileNotFoundError("Error in marker_reduction.py: AFLAP_tmp/LA.txt not found.")

    for G_info in get_LA_info():
        G, LO, UP, P0 = G_info

        with open(f"AFLAP_tmp/05/{G}_m{kmer}_L{LO}_U{UP}_{P0}.Filtered.Genotypes.MarkerID.tsv", 'r') as ftsv:
            MARKER_COUNT = len(ftsv.readlines())

        print(f"{MARKER_COUNT} markers detected in genotype file.")
        if (MARKER_COUNT > max_markers):
            print(f"Pruning til marker count is {max_markers}...")
            with open(f"AFLAP_tmp/05/{G}_m{kmer}_L{LO}_U{UP}_{P0}.Filtered.Genotypes.MarkerID.tsv", 'r+') as ftsv:
                lines = ftsv.readlines()
                lines = random.sample(lines, k=max_markers)
                ftsv.seek(0)
                ftsv.writelines(lines)
