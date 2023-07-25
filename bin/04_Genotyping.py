import argparse
import pandas as pd
import os
import subprocess

from get_LA_info import get_LA_info

#################################################
#	A Python script to call genotypes of progeny using markers derived from a parent and progeny JELLYFISH hashes.
#	For optimal calls, a k-mer should be observed twice.
#################################################

def create_count(count_file:str, prog:str, f_type:str, parent:str, kmer:int, lo:int, up:int, p0:str)->None:
    if os.path.exists(count_file) and os.path.getsize(count_file):
        print(f"\t\t\tCount for {prog} detected. Skipping")
        return

    jf_out = subprocess.run(f"jellyfish query -s AFLAP_tmp/03/F0Markers/{parent}_m{kmer}_MARKERS_L{lo}_U{up}_{p0}.fa AFLAP_tmp/01/{f_type}Count/{prog}.jf{kmer}",
                            shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
    with open(count_file, 'w') as f:
        for line in jf_out:
            if len(line): f.write(f"{line}\n")

    if not os.path.getsize(count_file):
        exit(f"An error occurred: Count file for {prog} was not created properly.")
    print(f"\t\t\tCount for {prog} created.")

def create_call(call_file:str, count_file:str, loci_seqs:set, prog:str, low_cov:int, f_type:str, sex:str)->None:
    if os.path.exists(call_file) and os.path.getsize(call_file):
        print(f"\t\t\tCall for {prog} detected. Skipping.")
        return

    with open(count_file, 'r') as fcount, open(call_file, 'w') as fcall:
        for line in fcount:
            line = line.strip().split()

            if f_type == "F2" and line[0] in loci_seqs: fcall.write("2\n")
            elif int(line[1]) >= low_cov: fcall.write("1\n")
            else: fcall.write("0\n")

    if not os.path.getsize(call_file):
        exit(f"An error occurred: Call file for {prog} was not created properly.")
    print(f"\t\t\tCall for {prog} created.")

def genotype_jfq(kmer:str, LowCov:str, G_info:tuple, f_type:str)->list:
    G, LO, UP, P0, SEX = G_info

    print(f"\tWorking on {f_type}...")

    # add progeny of parent to list
    ped_file = f"AFLAP_tmp/Pedigree_{f_type}.txt"
    if not os.path.exists(ped_file):
        exit(f"An error occurred: {ped_file} not found. Rerun AFLAP.py")

    prog_list = list()
    prog_df = pd.read_csv(f"AFLAP_tmp/Pedigree_{f_type}.txt", sep='\t')
    prog_list = prog_df[(prog_df["MP"].astype(str) == G) | (prog_df["FP"].astype(str) == G)]["Individual"].unique().tolist()
    if not len(prog_list):
        print(f"\t\tNo progeny of {G} found among given {f_type}.")

    # get loci sequences if analyzing 2nd generation type
    if f_type == "F2":
        same_loci_seqs = pd.read_csv("AFLAP_tmp/03/SimGroups/identical_loci.txt", sep='\t')
        set_of_loci_seqs = set(same_loci_seqs[f"{SEX.capitalize()} Sequence"].to_list())

    # perform jellyfish query
    for prog in prog_list:
        print(f"\t\tCreating Count and Call for {prog}...")
        if not os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{prog}.jf{kmer}"):
            exit(f"An error occurred: {prog} not detected among {f_type} progeny. Rerun 01_JELLYFISH.py.")

        count_file = f"AFLAP_tmp/04/Count/{prog}_{G}_m{kmer}_L{LO}_U{UP}_{P0}.txt"
        create_count(count_file, prog, f_type, G, kmer, LO, UP, P0)

        call_file = f"AFLAP_tmp/04/Call/{prog}_{G}_m{kmer}_L{LO}_U{UP}_{P0}.txt"
        create_call(call_file, count_file, set_of_loci_seqs, prog, int(LowCov), f_type, SEX)

    return prog_list

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
        G, LO, UP, P0, SEX = G_info

        # check if marker exists
        marker_file = f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa"
        if not os.path.exists(marker_file):
            exit(f"An error occurred: {marker_file} not found. Rerun 03_ObtainMarkers.py.")
        with open(marker_file, 'r') as fmark:
            m_count = 0
            for m in fmark:
                if m.startswith('>'):
                    m_count += 1
            print(f"\t{m_count} markers identified in {marker_file}. These will be surveyed against progeny.")

        # run jellyfish on F1 and F2 individuals who descend from F0
        prog_list = genotype_jfq(args.kmer, args.LowCov, G_info, "F1")
        prog_list += genotype_jfq(args.kmer, args.LowCov, G_info, "F2")

        # extract info from marker file
        marker_file = f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa"
        with open(marker_file, 'r') as f:
            head_list = list()
            seq_list = list()

            while True:
                head = f.readline().strip().replace('>', '')
                if not head: break

                seq = f.readline().strip()

                head_list.append(head)
                seq_list.append(seq)

            # check if head_list and seq_list are same size
            if len(head_list) != len(seq_list):
                exit(f"An error occurred: {marker_file} not extracted properly. Make sure every marker pairs with a sequence.")

        # get data
        data = {"MarkerSequence": seq_list, "MarkerID": head_list}
        for prog in prog_list:
            with open(f"AFLAP_tmp/04/Call/{prog}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt", 'r') as fcall:
                b_vals = list()
                for b_val in fcall: b_vals.append(b_val.strip())
            data[prog] = b_vals
        marker_df = pd.DataFrame(data=data)

        # split marker sequence and value and reorder
        marker_df[["MarkerID", "MarkerLength"]] = marker_df["MarkerID"].str.split('_', expand=True)
        marker_df = marker_df.reindex(columns=["MarkerSequence", "MarkerID", "MarkerLength"] + list(marker_df.columns[2:-1]))

        # create tsv file
        marker_df.to_csv(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t', index=False)

        if not os.path.exists(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv"):
            exit("An error occurred: Genotypes.MarkerID.tsv was not made.")
        print(f"\tGenotypes.MarkerID.tsv for {G} has been created.")
