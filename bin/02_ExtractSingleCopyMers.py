import argparse
import os
import subprocess

from get_LA_info import get_LA_info
from histoplot import histoplot

#################################################
#	A Python script to generate and plot histograms of parental hashes.
#	It will try to calculate peaks, if the user does not define them in the pedigree file, though this may be error prone.
#	Finally, it extracts k-mers it estimates to be single copy.
#################################################

def create_histogram(G:int, kmer:int)->None:
    if os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo") and os.path.getsize(f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo"):
        print(f"\tHistogram for {G} detected. Skipping.")
    else:
        if os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo") and not os.path.getsize(f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo"):
            print(f"\tEmpty histogram for {G} found. Deleting...")
            subprocess.run(f"rm AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo", shell=True)
            print(f"\tCreating new histogram for {G}...")
        else:
            print(f"\tCreating histogram for {G}...")

        subprocess.run(f"jellyfish histo AFLAP_tmp/01/F0Count/{G}.jf{kmer} > AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo",
                       shell=True, executable="/bin/bash")

        # check if jellyfish histo worked properly
        if not os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo"):
            raise FileNotFoundError(f"Jellyfish did not produce AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo.")
        if not os.path.getsize(f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo"):
            raise ValueError(f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo is empty.")

        print(f"\tHistogram for {G} generated.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ExtractingSingleCopyMers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()

    #  make directories
    os.makedirs("AFLAP_Results/Plots", exist_ok=True)
    os.makedirs("AFLAP_tmp/02/F0Histo", exist_ok=True)

    # get/make histograms
    try:
        print("Generating F0 histograms to undergo linkage analysis...")
        if not os.path.exists("AFLAP_tmp/LA.txt"):
            raise FileNotFoundError("AFLAP_tmp/LA.txt not found. Rerun 01_JELLYFISH.py.")
        list_of_Gs = get_LA_info()
        for G_info in list_of_Gs:
            G, LO, UP, P0 = G_info
            mer_count = 0

            create_histogram(G, args.kmer)

            print(f"\t\t{G} Bounds:\n" +
                f"\t\t\tLower: {LO}\n" +
                f"\t\t\tUpper: {UP}\n")

            # extract k-mers
            print(f"\tExtracting {args.kmer}-mers from {G}:")
            if os.path.exists(f"AFLAP_tmp/02/{G}_m{args.kmer}_L{LO}_U{UP}.fa"):
                print(f"\t\t{args.kmer}-mers for {G} detected. Skipping.")
                histo_same = True
            else:
                print(f"\t\tRunning jellyfish dump for {G}...")
                cmd = f"jellyfish dump -U {UP} -L {LO} -o AFLAP_tmp/02/{G}_m{args.kmer}_L{LO}_U{UP}.fa AFLAP_tmp/01/F0Count/{G}.jf{args.kmer}"
                subprocess.run(cmd, shell=True, executable="/bin/bash")
                histo_same = False

            # counting k-mers
            print(f"\tCounting number of {args.kmer}-mers for {G}...")
            print(f"\t\t{int(os.path.getsize(f'AFLAP_tmp/02/{G}_m{args.kmer}_L{LO}_U{UP}.fa') / 2)} {args.kmer}-mers counted for {G}")

            # create histo.png
            if os.path.exists(f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png") and histo_same:
                print(f"\tHistogram for {G} detected. Skipping.")
            else:
                if not histo_same:
                    print(f"\tMaking new histogram for {G}...")
                else:
                    print(f"\t Making histogram for {G}...")
                histoplot(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo", LO, UP, f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png")

                # check if histogram had been built
                if not os.path.exists(f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png"):
                    raise FileNotFoundError(f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png not found.")
                print(f"\tHistogram for {G} constructed.")

    except Exception as e:
        print(f"Error in 02_ExtractSingleCopyMers.py: {e}")
        exit(1)
