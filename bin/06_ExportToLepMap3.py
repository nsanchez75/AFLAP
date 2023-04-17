import argparse
import os

import get_LA_info as gli

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ExportToLepMap3', description="A script to export the genotype table to LepMap3.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()


    # 1. make directories
    os.makedirs("AFLAP_tmp/06", exist_ok=True)

    # 2. get parents to run LepMap3 on
    if (not os.path.exists("AFLAP_tmp/01/LA.txt")):
        raise FileNotFoundError("AFLAP_tmp/01/LA.txt not found. Rerun pipeline.")
    list_of_parents = gli.get_LA_info(f"AFLAP_tmp/01/LA.txt", f"AFLAP_tmp/01/Crosses.txt")
    for info in list_of_parents:
        G, LO, UP, P0 = info
        print(f"{G}\n{LO}\n{UP}\n{P0}")
    exit(0)

    with open("AFLAP_tmp/01/LA.txt", 'r') as fla:
        for p in fla:
            # extract info about analyzed parent
            p = p.strip().split()
            G = p[0]
            LO = p[1]
            UP = p[2]

            if (not os.path.exists(f"AFLAP_tmp/03/{G}_CrossedTo.txt")):
                raise FileNotFoundError(f"AFLAP_tmp/03/{G}_CrossedTo.txt not found. Rerun pipeline.")
            with open(f"AFLAP_tmp/03/{G}_CrossedTo.txt", 'r') as fct:
                p0 = []
                for op in fct:
                    p0.append(op.strip())
                p0 = '_'.join(p0)
            
            if (not os.path.exists(f"AFLAP_tmp/05/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.MarkerID.Filtered.tsv")):
                raise FileNotFoundError("Filtered .tsv file not found. Rerun 05_ObtainSegStats.py.")
            