import pandas as pd
import os
import subprocess

def jellyfish_count(kmer:str, threads:str, f_type:str)->None:
    if not os.path.getsize(f"AFLAP_tmp/Pedigree_{f_type}.txt"):
        print(f"No {f_type} individuals detected.")
        return
    print(f"Performing jellyfish {f_type} count:")

    ped_df = pd.read_csv(f"AFLAP_tmp/Pedigree_{f_type}.txt", sep='\t',
                         names=['Individual', 'Generation', 'File', 'LB', 'UB'])

    for ind in ped_df['Individual'].unique():
        # check if individual's count exists
        if os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
            if not os.path.getsize(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
                raise ValueError(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} is empty. " +
                                 f"Delete it and restart.")
            print(f"\tHash detected for {ind}. Skipping.")
            continue

        print(f"\tRunning jellyfish count for {ind}...")
        ind_df = ped_df.where(ped_df['Individual'] == ind).dropna()

        # get files to pass into jellyfish count
        jfin = list()
        for file in ind_df['File'].unique():
            if not os.path.exists(file):
                raise FileNotFoundError(f"{file} for {ind} not found.")
            jfin.append(file)
        jfin = ' '.join(jfin)

        # run jellyfish count
        cmd = f"jellyfish count -m {kmer} -C -s 1G -t {threads} -o AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} <(zcat {jfin})"
        subprocess.run(cmd, shell=True, executable="/bin/bash")

        # check if jellyfish worked
        if not os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
            raise FileNotFoundError(f"Jellyfish for {ind} did not complete.")
        print(f"\tHash for {ind} performed.")

    # with open(f"AFLAP_tmp/Pedigree_{f_type}.txt", 'r') as f:
    #     for ind in f:
    #         ind = ind.strip().split()[0]

    #         # check if individual's count exists
    #         if os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
    #             if not os.path.getsize(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
    #                 raise ValueError(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} is empty. " +
    #                                  f"Delete it and restart.")

    #             print(f"\tHash detected for {ind}. Skipping.")
    #         else:
    #             print(f"\tRunning jellyfish count for {ind}...")

    #             reads = []
    #             with open(f"AFLAP_tmp/Pedigree_{f_type}.txt", 'r') as fin:
    #                 for line in fin:
    #                     line = line.strip().split()
    #                     if line[0] == ind:
    #                         if not os.path.exists(line[2]):
    #                             raise FileNotFoundError(f"{line[2]} for {ind} not found.")
    #                         reads.append(line[2])

    #             # run jellyfish command
    #             jfin = ' '.join(reads)
    #             cmd = f"jellyfish count -m {kmer} -C -s 1G -t {threads} -o AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} <(zcat {jfin})"
    #             subprocess.run(cmd, shell=True, executable="/bin/bash")

    #             # check if jellyfish worked
    #             if not os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
    #                 raise FileNotFoundError(f"Jellyfish for {ind} did not complete.")
    #             else:
    #                 print(f"\tHash for {ind} performed.")

    print(f"Jellyfish {f_type} count complete.")