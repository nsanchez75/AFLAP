import os
import subprocess
import sys

def genotype_jfq(kmer:str, LOD:str, parent:str, lo:str, up:str, p0:str, f_type:str)->list:
    # declare which F type is being worked on
    print(f"\tWorking on {f_type}...")

    # define list variable
    h_list = []

    # add progeny of parent to list
    if os.path.exists(f"AFLAP_tmp/Pedigree_{f_type}.txt"):
        with open(f"AFLAP_tmp/Pedigree_{f_type}.txt") as f:
            prog_set = set()
            for prog in f:
                prog = prog.strip().split()

                if parent in {prog[3], prog[4]} and prog[0] not in prog_set:
                    h_list.append(prog[0])
                    prog_set.add(prog[0])

        # check if any progeny found
        if not len(h_list):
            print(f"Error in genotype_jf.py: No progeny of {parent} found among given {f_type}.")
            sys.exit(1)

        # perform jellyfish query if necessary
        for h in h_list:
            if not os.path.exists(f"AFLAP_tmp/01/{f_type}Count/{h}.jf{kmer}"):
                print(f"Error in genotype_jf.py: {h} not detected among {f_type} progeny. Rerun 01_JELLYFISH.py.")
                sys.exit(1)

            # create counts of progeny
            if os.path.exists(f"AFLAP_tmp/04/Count/{h}_{parent}_m{kmer}_L{lo}_U{up}_{p0}.txt"):
                print(f"\tCount for {h} detected. Skipping")
            else:
                jf_cmd = f"jellyfish query -s AFLAP_tmp/03/F0Markers/{parent}_m{kmer}_MARKERS_L{lo}_U{up}_{p0}.fa AFLAP_tmp/01/{f_type}Count/{h}.jf{kmer}"
                jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                with open(f"AFLAP_tmp/04/Count/{h}_{parent}_m{kmer}_L{lo}_U{up}_{p0}.txt", 'w') as f:
                    for line in jf_out:
                        # disregard empty lines FIXME: determine why this happens
                        if not len(line): continue

                        f.write(f"{line}\n")

            # create calls of progeny
            if os.path.exists(f"AFLAP_tmp/04/Call/{h}_{parent}_m{kmer}_L{lo}_U{up}_{p0}.txt"):
                print(f"\tCall for {h} detected. Skipping.")
            else:
                with open(f"AFLAP_tmp/04/Count/{h}_{parent}_m{kmer}_L{lo}_U{up}_{p0}.txt", 'r') as fcount, open(f"AFLAP_tmp/04/Call/{h}_{parent}_m{kmer}_L{lo}_U{up}_{p0}.txt", 'w') as fcall:
                    for line in fcount:
                        line = line.strip().split()

                        if int(line[1]) >= int(LOD): fcall.write(f"1\n")
                        else: fcall.write(f"0\n")

    return h_list