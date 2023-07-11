import pandas as pd
import os
import subprocess

def jellyfish_count(kmer:str, threads:str, f_type:str)->None:
    if not os.path.getsize(f"AFLAP_tmp/Pedigree_{f_type}.txt"):
        print(f"No {f_type} individuals detected.")
        return
    print(f"Performing jellyfish {f_type} count:")

    print("current directory:", os.path.curdir)

    ped_df = pd.read_csv(f"AFLAP_tmp/Pedigree_{f_type}.txt", sep='\t')

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
            if not os.path.exists(file):
                exit(f"An error occurred: {file} for {ind} not found.")
            jfin.append(file)
        jfin = ' '.join(jfin)

        subprocess.run(args=f"jellyfish count -m {kmer} -C -s 1G -t {threads} -o AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} <(zcat {jfin})",
                       shell=True, executable="/bin/bash")
        # check if jellyfish worked
        if not os.path.exists(jc_file) or not os.path.getsize(jc_file):
            exit(f"An error occurred: Jellyfish for {ind} did not complete.")
        print(f"\tHash for {ind} performed.")

    print(f"Jellyfish {f_type} count complete.")