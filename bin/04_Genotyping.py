import argparse
import glob
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
        [f.write(f"{line}\n") if (len(line)) else None for line in jf_out]

    if not os.path.getsize(count_file):
        exit(f"An error occurred: Count file for {prog} was not created properly.")
    print(f"\t\t\tCount for {prog} created.")

def create_call(call_file:str, count_file:str, prog:str, low_cov:int, f_type:str, sex:str)->None:
    if os.path.exists(call_file) and os.path.getsize(call_file):
        print(f"\t\t\tCall for {prog} detected. Skipping.")
        return

    with open(count_file, 'r') as fcount, open(call_file, 'w') as fcall:
        for line in fcount:
            line = line.strip().split()
            if int(line[1]) >= low_cov: fcall.write("1\n")
            else:                       fcall.write("0\n")

    if not os.path.getsize(call_file):
        exit(f"An error occurred: Call file for {prog} was not created properly.")
    print(f"\t\t\tCall for {prog} created.")

def genotype_jfq(kmer:str, LowCov:str, G_info:tuple, f_type:str)->list:
    ped_file = f"AFLAP_tmp/Pedigree_{f_type}.txt"
    if not os.path.exists(ped_file):
        exit(f"An error occurred: {ped_file} not found. Rerun AFLAP.py")
    if not os.path.getsize(ped_file):
        print(f"No {f_type} progeny detected. Skipping.")
        return

    print(f"\tWorking on {f_type}...")
    G, LO, UP, P0, SEX = G_info

    # add progeny of parent to list
    prog_list = list()
    prog_df = pd.read_csv(f"AFLAP_tmp/Pedigree_{f_type}.txt", sep='\t')
    prog_list = prog_df[(prog_df["MP"].astype(str) == G) | (prog_df["FP"].astype(str) == G)]["Individual"].unique().tolist()
    if not len(prog_list):
        print(f"\t\tNo progeny of {G} found among given {f_type}.")

    # perform jellyfish query
    for prog in prog_list:
        print(f"\t\tCreating Count and Call for {prog}...")
        if not os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{prog}.jf{kmer}"):
            exit(f"An error occurred: {prog} not detected among {f_type} progeny. Rerun 01_JELLYFISH.py.")

        count_file = f"AFLAP_tmp/04/{f_type}/Count/{prog}_{G}_m{kmer}_L{LO}_U{UP}_{P0}.txt"
        create_count(count_file, prog, f_type, G, kmer, LO, UP, P0)
        call_file = f"AFLAP_tmp/04/{f_type}/Call/{prog}_{G}_m{kmer}_L{LO}_U{UP}_{P0}.txt"
        create_call(call_file, count_file, prog, int(LowCov), f_type, SEX)

    num_count_files = len(list(filter(lambda x: True if (x[22:].split('_')[1] == G) else False, glob.glob(f"AFLAP_tmp/04/{f_type}/Count/*.txt"))))
    num_call_files = len(list(filter(lambda x: True if (x[21:].split('_')[1] == G) else False, glob.glob(f"AFLAP_tmp/04/{f_type}/Call/*.txt"))))
    if not num_count_files == num_call_files:
        exit(f"An error occurred: Non-equivalent amount of Count and Call files were made for {f_type} progeny of {G}. Delete 04 directory in AFLAP_tmp and rerun 04_Genotyping.py.")

    return prog_list

def make_f2_genotype_table(marker_df:pd.DataFrame, ident_loci_df:pd.DataFrame, sex:str)->pd.DataFrame:
    # edit identical loci dataframe
    id_col_name = f"{sex.capitalize()} Sequence ID"
    specific_loci_df = ident_loci_df[[id_col_name, "Locus Sequence", "Locus Sequence ID"]]
    # TODO: fix warning
    specific_loci_df[id_col_name] = specific_loci_df[id_col_name].apply(lambda x: x.split('_')[0])

    # merge dataframes and replace marker ID with locus ID
    marker_df = marker_df.merge(specific_loci_df, left_on="MarkerID", right_on=id_col_name, how='inner')
    marker_df["MarkerSequence"] = marker_df["Locus Sequence"]
    marker_df["MarkerID"] = marker_df["Locus Sequence ID"]

    if sex == "male": marker_df.iloc[:, 3:] = marker_df.iloc[:, 3:].replace(['0', '1'], ['X', 'A'])
    else:             marker_df.iloc[:, 3:] = marker_df.iloc[:, 3:].replace(['0', '1'], ['X', 'B'])

    return marker_df.drop(columns=[id_col_name, "Locus Sequence", "Locus Sequence ID"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Genotyping', description="A script to genotype progeny")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-x', '--LowCov', type=int, default=2, help='Run with low coverage parameters.')
    args = parser.parse_args()

    list_of_Gs = get_LA_info()
    for f_type in ["F1", "F2"]:
        if not len(os.listdir(f"AFLAP_tmp/01/{f_type}Count")):
            print(f"No {f_type} progeny found. Skipping.")
            continue

        # make directories
        os.makedirs(f"AFLAP_tmp/04/{f_type}/Count", exist_ok=True)
        os.makedirs(f"AFLAP_tmp/04/{f_type}/Call", exist_ok=True)

        # check for markers
        for G_info in list_of_Gs:
            G, LO, UP, P0, SEX = G_info

            # check if marker exists
            marker_file = f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa"
            if not os.path.exists(marker_file) or not os.path.getsize(marker_file):
                exit(f"An error occurred: {marker_file} not valid. Rerun 03_ObtainMarkers.py.")
            print(f"\t{os.path.getsize(marker_file) // 2} markers identified in {marker_file}. These will be surveyed against progeny.")

            # get data for progeny
            print(f"\tGetting {f_type} progeny data for {G}...")
            prog_list = genotype_jfq(args.kmer, args.LowCov, G_info, f_type)
            if not prog_list:
                exit(f"An error occurred: Creating call and count of {f_type} progeny did not work.")

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

                if len(head_list) != len(seq_list):
                    exit(f"An error occurred: {marker_file} not extracted properly. Make sure every marker pairs with a sequence.")

            # create genotype table
            print(f"Creating a genotype table for {f_type} progeny of {G}...")
            data = {"MarkerSequence": seq_list, "MarkerID": head_list}
            for prog in prog_list:
                with open(f"AFLAP_tmp/04/{f_type}/Call/{prog}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt", 'r') as fcall:
                    b_vals = list()
                    for b_val in fcall: b_vals.append(b_val.strip())
                data[prog] = b_vals
            marker_df = pd.DataFrame(data=data)

            ## split marker sequence and value and reorder
            marker_df[["MarkerID", "MarkerLength"]] = marker_df["MarkerID"].str.split('_', expand=True)
            marker_df = marker_df.reindex(columns=["MarkerSequence", "MarkerID", "MarkerLength"] + list(marker_df.columns[2:-1]))

            if f_type == "F2":
                identical_loci_id_df = pd.read_csv("AFLAP_Results/IdenticalLoci.txt", sep='\t')
                if identical_loci_id_df.empty:
                    exit("An error occurred: AFLAP_Results/IdenticalLoci.txt not found. Rerun 03_ObtainMarkers.py.")

                marker_df = make_f2_genotype_table(marker_df, identical_loci_id_df, SEX)

            ## create tsv file
            marker_df.to_csv(f"AFLAP_tmp/04/{G}_{f_type}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t', index=False)

            if not os.path.exists(f"AFLAP_tmp/04/{G}_{f_type}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv"):
                exit("An error occurred: Genotypes.MarkerID.tsv was not made.")
            print(f"\t{f_type} Genotypes.MarkerID.tsv for {G} has been created.")
