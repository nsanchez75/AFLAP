import argparse
import os
import sys

import jellyfish_count as jc

###########################################################
#	Python script to be ran as part of the AFLAP pipeline.
#	This script will sanity check the pedigree file then run JELLYFISH on parents and progeny specified in the pedigree file.
#	The output will be stored in AFLAP_tmp.
#	The script will detect previous results and used them when able.
###########################################################

def main()->None:
    parser = argparse.ArgumentParser(prog='JELLYFISH', description="A script to obtain JELLYFISH hashes to be used when running AFLAP.")
    parser.add_argument('-m', '--kmer', default=31, help="K-mer size (optional). Default [31].")
    parser.add_argument('-t', '--threads', default=4, help="Threads for JELLYFISH counting (optional). Default [4].")
    args = parser.parse_args()


    # 1. indentify which parents will have a map constructed for them
    ## initialize variables
    parents    = set()
    la_count   = set()
    nola_count = set()

    ## perform LA analysis
    with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as fin, open("AFLAP_tmp/01/LA.txt", 'w') as fla, open("AFLAP_tmp/01/noLA.txt", 'w') as fnola:
        for line in fin:
            line = line.strip().split()
            if line[0] not in parents:
                if "NA" in [line[3], line[4]]:
                    fnola.write(f"{line[0]}\n")
                    parents.add(line[0])
                elif [isinstance(x, int) for x in [line[3], line[4]]]:
                    if line[3] > line[4]:
                        print("Error: cannot have a lower bound higher than an upper bound.")
                        sys.exit(1)

                    # write {parent} {lower bound} {upper bound} to LA.txt
                    fla.write(f"{line[0]} {line[3]} {line[4]}")
                    parents.add(line[0])
                else:
                    print("Error: invalid bound entry in Pedigree_F0.txt.")
                    sys.exit(1)


    # 2. perform k-mer counting
    ## make directories
    os.makedirs("AFLAP_tmp/01/F0_Count", exist_ok=True)
    os.makedirs("AFLAP_tmp/01/F1_Count", exist_ok=True)
    os.makedirs("AFLAP_tmp/01/F2_Count", exist_ok=True)

    ## parents
    jc.jellyfish_count(args.kmer, args.threads, "F0")
    ## F1 progeny
    jc.jellyfish_count(args.kmer, args.threads, "F1")
    ## F2 progeny
    # TODO: uncomment this when working on F2
    # jc.jellyfish_count(args.kmer, args.threads, "F2")


if __name__ == "__main__":
    main()