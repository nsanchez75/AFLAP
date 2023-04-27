import argparse
import os
import subprocess

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
    list_of_Gs = gli.get_LA_info()
    for G_info in list_of_Gs:
        G, LO, UP, P0 = G_info

        # check if tsv file for G exists
        if not os.path.exists(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv"):
            raise FileNotFoundError(f"Error: LepMap-ready genotype table fo {G} not found.")

        # begin LepMap3 analysis
        print(f"Initiating LepMap3 analysis on {G}...")
        ## skip if G had already been analyzed
        if os.path.exists(f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.txt"):
            print(f"\tPrevious results for {G} detected. Skipping.")
        else:
            # run LepMap3 - SeparateChromosomes2
            sc2_stdout = open(f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.txt", 'w')
            sc2_stderr = open(f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.stderr", 'w')
            subprocess.run(f"java -cp $CONDA_PREFIX/bin/lepmap3/ SeparateChromosomes2 lodLimit={args.LOD} numThreads={args.threads} data=AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv",
                           stdout=sc2_stdout, stderr=sc2_stderr)
            sc2_stdout.close(), sc2_stderr.close()

        print("continue coding 07?")