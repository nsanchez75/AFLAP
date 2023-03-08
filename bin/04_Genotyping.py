import argparse
import os
import pandas as pd
import subprocess

import genotype_jf as gj

def main()->None:
    parser = argparse.ArgumentParser(prog='Genotyping', description="A script to genotype progeny")
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()


    # 1. make directories
    os.makedirs("AFLAP_tmp/04/Count", exist_ok=True)
    os.makedirs("AFLAP_tmp/04/Call", exist_ok=True)


    # 2. check for markers
    with open("AFLAP_tmp/01/LA.txt", 'r') as fla:
        for line in fla:
            line = line.strip().split()

            # declare variables from LA.txt
            G  = line[0]    # F0 being analyzed
            LO = line[1]    # F0's lower bound
            UP = line[2]    # F1's upper bound

            print(f"Calling GT for {G}'s derived markers...")

            # identify file based on intermediate results
            with open(f"AFLAP_tmp/03/{G}_CrossedTo.txt", 'r') as fct:
                p0 = []
                for op in fct:
                    p0.append(op.strip())
                p0 = 'x'.join(p0)

            # check if marker exists
            if os.path.exists(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa"):
                with open(f"AFLAP_tmp/03/ParentalMarkers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa", 'r') as fmark:
                    m_count = 0
                    for m in fmark:
                        if m.startswith('>'):
                            m_count += 1

                    print(f"\t{m_count} markers identified in AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa. These will be surveyed against progeny.")
            else:
                print(f"Error in 04_Genotyping.py: AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa not found. Rerun 03_ObtainMarkers.py.")
                exit(1)

            # run jellyfish on F1 and F2 individuals who descend from F0
            gj.genotype_jfq(args.kmer, args.LOD, G, LO, UP, p0, "F1")
            gj.genotype_jfq(args.kmer, args.LOD, G, LO, UP, p0, "F2")

            # put contents into .tsv file
            with open(f"AFLAP_tmp/04/Count/{h}_{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.txt", 'r') as fcount, open(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.tsv", 'w') as fout:
                for line in fcount:
                    line = line.strip().split()

                    fcall = f"AFLAP_tmp/04/Call/*_{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.txt"
                    paste_output = os.popen(f"paste - {fcall}")
                    df = pd.DataFrame([i.strip().split() for i in paste_output.strip('\n') if i])


if __name__ == "__main__":
    main()