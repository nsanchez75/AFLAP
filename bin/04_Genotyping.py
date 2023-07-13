import argparse
import pandas as pd
import os
import subprocess

from get_LA_info import get_LA_info

#################################################
#	A Python script to call genotypes of progeny using markers derived from a parent and progeny JELLYFISH hashes.
#	For optimal calls, a k-mer should be observed twice.
#################################################

def genotype_jfq(kmer:str, LowCov:str, parent:str, lo:str, up:str, p0:str, f_type:str)->list:
    # declare which F type is being worked on
    print(f"\tWorking on {f_type}...")

    # add progeny of parent to list
    ped_file = f"AFLAP_tmp/Pedigree_{f_type}.txt"
    if not os.path.exists(ped_file):
        exit(f"An error occurred: {ped_file} not found. Rerun AFLAP.py")

    ped_df = pd.read_csv(ped_file, sep='\t')
    # TODO: if this works then refactor all '.loc' stuff
    print(ped_df)
    print(parent)
    prog_df = ped_df[(ped_df["MP"] == parent) | (ped_df["FP"] == parent)]
    print(prog_df)
    h_list = prog_df["Individual"].unique().tolist()

    # check if any progeny found
    if not len(h_list):
        print(f"\t\tNo progeny of {parent} found among given {f_type}.")

    # perform jellyfish query if necessary
    for h in h_list:
        print(f"\t\tCreating Count and Call for {h}...")
        if not os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{h}.jf{kmer}"):
            exit(f"An error occurred: {h} not detected among {f_type} progeny. Rerun 01_JELLYFISH.py.")

        # create counts of progeny
        count_file = f"AFLAP_tmp/04/Count/{h}_{parent}_m{kmer}_L{lo}_U{up}_{p0}.txt"
        if os.path.exists(count_file) and os.path.getsize(count_file):
            print(f"\t\t\tCount for {h} detected. Skipping")
        else:
            jf_cmd = f"jellyfish query -s AFLAP_tmp/03/F0Markers/{parent}_m{kmer}_MARKERS_L{lo}_U{up}_{p0}.fa AFLAP_tmp/01/{f_type}Count/{h}.jf{kmer}"
            jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
            with open(count_file, 'w') as f:
                for line in jf_out:
                    if not len(line): continue
                    f.write(f"{line}\n")
            print(f"\t\t\tCount for {h} created.")

        # create calls of progeny
        call_file = f"AFLAP_tmp/04/Call/{h}_{parent}_m{kmer}_L{lo}_U{up}_{p0}.txt"
        if os.path.exists(call_file) and os.path.getsize(call_file):
            print(f"\t\t\tCall for {h} detected. Skipping.")
        else:
            with open(count_file, 'r') as fcount, open(call_file, 'w') as fcall:
                for line in fcount:
                    line = line.strip().split()

                    if int(line[1]) >= int(LowCov): fcall.write(f"1\n")
                    else: fcall.write(f"0\n")
            print(f"\t\t\tCall for {h} created.")

    return h_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Genotyping', description="A script to genotype progeny")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-x', '--LowCov', type=int, default=2, help='Run with low coverage parameters.')
    args = parser.parse_args()

    # make directories
    os.makedirs("AFLAP_tmp/04/Count", exist_ok=True)
    os.makedirs("AFLAP_tmp/04/Call", exist_ok=True)

    # check for markers
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
            (f"An error occurred: AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa not found. Rerun 03_ObtainMarkers.py.")

        # initialize list consisting of F0's progeny
        h_list = []

        # run jellyfish on F1 and F2 individuals who descend from F0
        h_list += genotype_jfq(args.kmer, args.LowCov, G, LO, UP, P0, "F1")
        h_list += genotype_jfq(args.kmer, args.LowCov, G, LO, UP, P0, "F2")

        # extract info from MARKERS file
        with open(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa", 'r') as f:
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
                exit(f"An error occurred: AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa not extracted properly.")

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
            exit("An error occurred: Genotypes.MarkerID.tsv was not made.")
        else:
            print(f"\tGenotypes.MarkerID.tsv for {G} has been created.")
