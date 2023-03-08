import os
import subprocess
import sys

def jellyfish_count(kmer:str, threads:str, f_type:str)->None:
    print(f"Performing jellyfish {f_type} count:")

    with open(f"AFLAP_tmp/01/{f_type}.txt", 'r') as f:
        for ind in f:
            ind = ind.strip()

            # check if individual's count exists
            if os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
                if not os.path.getsize(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
                    print(f"Error: AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} is empty.")
                    sys.exit(1)

                print(f"\tHash detected for {ind}. Skipping.")
            else:
                print(f"\tRunning jellyfish count for {ind}")

                reads = []
                with open(f"AFLAP_tmp/Pedigree_{f_type}.txt", 'r') as fin:
                    for line in fin:
                        line = line.strip().split()
                        if line[0] == ind:
                            if not os.path.exists(line[1]):
                                print(f"Error: {line[1]} not found.")
                                sys.exit(1)
                            reads.append(line[2])

                # run jellyfish command
                jfin = ' '.join(reads)
                cmd = f"jellyfish count -m {kmer} -C -s 1G -t {threads} -o AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer} <(zcat {jfin})"
                subprocess.run(cmd, shell=True, executable="/bin/bash")

                # check if jellyfish worked
                if not os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{ind}.jf{kmer}"):
                    print(f"Error: Jellyfish for {ind} did not complete.")
                    sys.exit(1)
                else:
                    print(f"\tHash for {ind} performed.")
