import argparse
import pandas as pd
import os

from get_LA_info import get_LA_info
from genotype_jf import genotype_jfq

#################################################
#	A Python script to call genotypes of progeny using makrers derived from a parent and progeny JELLYFISH hashes.
#	For optimal calls, a k-mer should be observed twice.
#################################################

def tsv_sort(line:str)->str:
    line_fields = line.strip().split()
    return str(line_fields[0])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Genotyping', description="A script to genotype progeny")
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()


    # make directories
    os.makedirs("AFLAP_tmp/04/Count", exist_ok=True)
    os.makedirs("AFLAP_tmp/04/Call", exist_ok=True)


    # check for markers
    try:
        list_of_Gs = get_LA_info()
        for G_info in list_of_Gs:
            G, LO, UP, P0 = G_info

            # check if marker exists
            if os.path.exists(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa"):
                with open(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa", 'r') as fmark:
                    m_count = 0
                    for m in fmark:
                        if m.startswith('>'):
                            m_count += 1

                    print(f"\t{m_count} markers identified in AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa. These will be surveyed against progeny.")
            else:
                raise FileNotFoundError(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa not found. Rerun 03_ObtainMarkers.py.")

            # initialize list consisting of F0's progeny
            h_list = []

            # run jellyfish on F1 and F2 individuals who descend from F0
            h_list += genotype_jfq(args.kmer, args.LOD, G, LO, UP, P0, "F1")
            h_list += genotype_jfq(args.kmer, args.LOD, G, LO, UP, P0, "F2")

            # extract info from MARKERS file
            with open(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa") as f:
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
                    raise ValueError(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa not extracted properly.")

            # get data
            data = {"MarkerSequence": seq_list, "MarkerID": head_list}
            for h in h_list:
                with open(f"AFLAP_tmp/04/Call/{h}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt", 'r') as fcall:
                    b_vals = []
                    for b_val in fcall: b_vals.append(b_val.strip())
                data[h] = b_vals

            matrix = pd.DataFrame(data=data)

            # split marker sequence and value and reorder
            matrix[["MarkerID", "MarkerLength"]] = matrix["MarkerID"].str.split('_', expand=True)
            matrix = matrix.reindex(columns=["MarkerSequence", "MarkerID", "MarkerLength"] + list(matrix.columns[2:-1]))

            # create tsv file
            matrix.to_csv(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t', index=False)

            # check tsv file status
            if not os.path.exists(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv"):
                raise FileNotFoundError("Genotypes.MarkerID.tsv was not made.")
            else:
                print(f"\tGenotypes.MarkerID.tsv for {G} has been created.")

    except Exception as e:
        print(f"Error in 04_Genotyping.py: {e}")
        exit(1)