import argparse
import pandas as pd
import os
import sys

import genotype_jf as gj

def tsv_sort(line:str)->str:
    line_fields = line.strip().split()
    return str(line_fields[0])

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
                p0 = '_'.join(p0)

            # check if marker exists
            if os.path.exists(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa"):
                with open(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa", 'r') as fmark:
                    m_count = 0
                    for m in fmark:
                        if m.startswith('>'):
                            m_count += 1

                    print(f"\t{m_count} markers identified in AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa. These will be surveyed against progeny.")
            else:
                print(f"Error in 04_Genotyping.py: AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa not found. Rerun 03_ObtainMarkers.py.")
                exit(1)

            # initialize list consisting of F0's progeny
            h_list = []

            # run jellyfish on F1 and F2 individuals who descend from F0
            h_list += gj.genotype_jfq(args.kmer, args.LOD, G, LO, UP, p0, "F1")
            h_list += gj.genotype_jfq(args.kmer, args.LOD, G, LO, UP, p0, "F2")

            # extract info from MARKERS file
            with open(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa") as f:
                head_list = []
                seq_list = []

                while True:
                    head = f.readline().strip().replace('>', '')
                    if not head: break

                    seq = f.readline().strip()

                    head_list.append(head)
                    seq_list.append(seq)

                # check if head_list and seq_list are same size
                if len(head_list) != len(seq_list):
                    print(f"Error in 04_Genotyping.py: AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa not extracted properly.")
                    sys.exit(1)

            # get data
            data = {"MarkerID": seq_list, "MarkerSequence": head_list}
            for h in h_list:
                with open(f"AFLAP_tmp/04/Call/{h}_{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.txt", 'r') as fcall:
                    b_vals = []
                    for b_val in fcall: b_vals.append(b_val.strip())
                data[h] = b_vals

            matrix = pd.DataFrame(data=data)

            # create tsv file
            matrix.to_csv(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.MarkerID.tsv", sep='\t')
            
            # check tsv file status
            if not os.path.exists(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.MarkerID.tsv"):
                print("Error in 04_Genotyping.py: Genotypes.MarkerID.tsv was not made.")
                sys.exit(1)
            else:
                print(f"Genotypes.MarkerID.tsv for {G} has been created.")


if __name__ == "__main__":
    main()