import argparse
import os

import get_LA_info as gli

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='LepMap3', description='A script to run LepMap3 and produce a genetic map which can be aligned to a genome assembly.')
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for JELLYFISH counting (optional). Default [4].')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()

    # check if necessary files exist
    if not os.path.exists("AFLAP_tmp/01/LA.txt"):
        raise FileNotFoundError("Error: AFLAP_tmp/01/LA.txt not found.")
    
    # make directory
    os.makedirs(f"AFLAP_Results/LOD{args.LOD}", exist_ok=True)
    
    # run LepMap3 on parents
    gli.get_LA_info("AFLAP_tmp/")