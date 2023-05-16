import argparse
import os

import jellyfish_count as jc

###########################################################
#	Python script to be ran as part of the AFLAP pipeline.
#	This script will sanity check the pedigree file then run JELLYFISH on parents and progeny specified in the pedigree file.
#	The output will be stored in AFLAP_tmp.
#	The script will detect previous results and used them when able.
###########################################################

if __name__ == "__main__":
    try:
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
        jc.jellyfish_count(args.kmer, args.threads, "F0")
        ## F1 progeny
        jc.jellyfish_count(args.kmer, args.threads, "F1")
        ## F2 progeny
        # TODO: uncomment this when working on F2
        jc.jellyfish_count(args.kmer, args.threads, "F2")

    except Exception as e:
        print(f"Error when running 01_JELLYFISH {e}")
        exit(1)