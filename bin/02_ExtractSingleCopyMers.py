import argparse
import os
import subprocess
import sys

from histoplot import histoplot

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ExtractingSingleCopyMers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes.")
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()


    # 1. make directories
    os.makedirs("AFLAP_Results/Plots", exist_ok=True)
    os.makedirs("AFLAP_tmp/02/F0Histo", exist_ok=True)


    # 2. get/make histograms
    print("Generating F0 histograms to undergo linkage analysis...")
    if not os.path.exists("AFLAP_tmp/01/LA.txt"):
        print("Error in 02_ExtractSingleCopyMers.py: AFLAP_tmp/01/LA.txt not found. Rerun 01_JELLYFISH.py.")
    with open("AFLAP_tmp/01/LA.txt", 'r') as fla:
        for p in fla:
            p = p.strip().split()

            # initialize variables
            G = int(p[0])
            LO = int(p[1])
            UP = int(p[2])
            mer_count = 0

            if os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo") and os.path.getsize(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo"):
                print(f"\tHistogram for {G} detected. Skipping.")
            else:
                if os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo") and not os.path.getsize(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo"):
                    print(f"\tEmpty histogram for {G} found. Deleting...")
                    subprocess.run(f"rm AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo", shell=True)
                    print(f"\tCreating new histogram for {G}...")
                else:
                    print(f"\tCreating histogram for {G}...")
                
                cmd = f"jellyfish histo AFLAP_tmp/01/F0Count/{G}.jf{args.kmer} > AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo"
                subprocess.run(cmd, shell=True, executable="/bin/bash")

                if not os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo"):
                    print(f"Error in 02_ExtractSingleCopyMers.py: jellyfish did not produce AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo.")
                    sys.exit(1)
                else:
                    # check if jellyfish histo worked properly
                    with open(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo", 'r') as f:
                        if len(f.readlines()) == 0:
                            print(f"Error in 02_ExtractSingleCopyMers.py: AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo is empty.")
                            sys.exit(1)

                    print(f"\tHistogram for {G} generated.")

            print(f"\t\t{G} Bounds:\n" +
                  f"\t\t\tLower: {LO}\n" +
                  f"\t\t\tUpper: {UP}")

            # extract k-mers
            print()
            print(f"\tExtracting {args.kmer}-mers from {G}:")
            if os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa"):
                print(f"\t\t{args.kmer}-mers for {G} detected. Skipping.")
                histo_same = True
            else:
                print(f"\t\tRunning jellyfish dump for {G}...")
                cmd = f"jellyfish dump -U {UP} -L {LO} -o AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa AFLAP_tmp/01/F0Count/{G}.jf{args.kmer}"
                subprocess.run(cmd, shell=True, executable="/bin/bash")
                histo_same = False

            # counting k-mers
            print(f"\tCounting number of {args.kmer}-mers for {G}...")
            with open(f"AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa", 'r') as fkmers:
                for line in fkmers:
                    if line.startswith('>'): mer_count += 1
                print(f"\t\t{mer_count} {args.kmer}-mers counted for {G}.")

            # create histo.png
            if os.path.exists(f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png") and histo_same:
                print(f"\t\tHistogram for {G} detected. Skipping.")
            else:
                if not histo_same:
                    print(f"\t\tMaking new histogram for {G}...")
                else:
                    print(f"\t\t Making histogram for {G}...")
                histoplot(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo", LO, UP, f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png")
                # cmd = f"Rscript bin/HistoPlot.R AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo {LO} {UP} AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png"
                # subprocess.run(cmd, shell=True)

                # check if histogram had been built
                if not os.path.exists(f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png"):
                    print(f"Error in 02_ExtractSingleCopyMers.py: AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png not found.")
                    sys.exit(1)
                else:
                    print(f"\t\tHistogram for {G} constructed.")
