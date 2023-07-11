import argparse
import os
import pandas as pd
import subprocess

###########################################################
#	A Python script to be ran as part of the AFLAP pipeline.
#	This script will sanity check the pedigree file then run JELLYFISH on parents and progeny specified in the pedigree file.
#	The output will be stored in AFLAP_tmp.
#	The script will detect previous results and used them when able.
###########################################################

def jellyfish_count(kmer:str, threads:str, f_type:str, ped_df:pd.DataFrame)->None:
    if not os.path.getsize(f"AFLAP_tmp/Pedigree_{f_type}.txt"):
        print(f"No {f_type} individuals detected.")
        return
    print(f"Performing jellyfish {f_type} count:")

    # run jellyfish count on each individual
    for ind in ped_df['Individual'].unique():
        jc_file = f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"
        if os.path.exists(jc_file) and os.path.getsize(jc_file):
            print(f"\tHash detected for {ind}. Skipping.")
            continue
        print(f"\tRunning jellyfish count for {ind}...")
        ind_df = ped_df[ped_df["Individual"] == ind].dropna()

        # get files to pass into jellyfish count
        jfin = list()
        for file in ind_df["Path"].unique():
            file = os.getcwd() + '/' + file
            print(file)
            # print("current path: " + os.getcwd()) # FIXME: get current path
            # print("the path " + file + ' exists? -> ' + str(os.path.exists(file)))
            if not os.path.exists(file):
                exit(f"An error occurred: {file} for {ind} not found. Make sure that your files' root path is in your current directory.")
            jfin.append(file)
        jfin = ' '.join(jfin)

        subprocess.run(args=f"jellyfish count -m {kmer} -C -s 1G -t {threads} -o AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} <(zcat {jfin})",
                       shell=True, executable="/bin/bash")
        # check if jellyfish worked
        if not os.path.exists(jc_file) or not os.path.getsize(jc_file):
            exit(f"An error occurred: Jellyfish for {ind} did not complete.")
        print(f"\tHash for {ind} performed.")

    print(f"Jellyfish {f_type} count complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='JELLYFISH', description="A script to obtain JELLYFISH hashes to be used when running AFLAP.")
    parser.add_argument('-m', '--kmer', default=31, help="K-mer size (optional). Default [31].")
    parser.add_argument('-t', '--threads', default=4, help="Threads for JELLYFISH counting (optional). Default [4].")
    args = parser.parse_args()

    # make directories
    os.makedirs("AFLAP_tmp/01/F0Count", exist_ok=True)
    os.makedirs("AFLAP_tmp/01/F1Count", exist_ok=True)
    os.makedirs("AFLAP_tmp/01/F2Count", exist_ok=True)

    # perform jellyfish counting
    ## parents
    ped_df = pd.read_csv(f"AFLAP_tmp/Pedigree_F0.txt", sep='\t')
    jellyfish_count(args.kmer, args.threads, "F0", ped_df)
    ## F1 progeny
    ped_df = pd.read_csv(f"AFLAP_tmp/Pedigree_F1.txt", sep='\t')
    jellyfish_count(args.kmer, args.threads, "F1", ped_df)
    ## F2 progeny
    ped_df = pd.read_csv(f"AFLAP_tmp/Pedigree_F2.txt", sep='\t')
    jellyfish_count(args.kmer, args.threads, "F2", ped_df)
