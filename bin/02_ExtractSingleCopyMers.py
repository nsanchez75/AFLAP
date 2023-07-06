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
    histo_file = f"AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo"
    if os.path.exists(histo_file) and os.path.getsize(histo_file):
        print(f"\tHistogram for {G} detected. Skipping.")
        return
    print(f"\tCreating histogram for {G}...")

    subprocess.run(f"jellyfish histo AFLAP_tmp/01/F0Count/{G}.jf{kmer} > AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo",
                    shell=True, executable="/bin/bash")
    # check if jellyfish histo worked properly
    if not os.path.exists(histo_file) or not os.path.getsize(histo_file):
        exit(f"An error occurred: Jellyfish did not create AFLAP_tmp/02/F0Histo/{G}.{kmer}.histo properly.")
    print(f"\tHistogram for {G} generated.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ExtractingSingleCopyMers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()

    #  make directories
    os.makedirs("AFLAP_Results/Plots", exist_ok=True)
    os.makedirs("AFLAP_tmp/02/F0Histo", exist_ok=True)

    # get/make histograms
    print("Generating F0 histograms to undergo linkage analysis...")
    if not os.path.exists("AFLAP_tmp/LA.txt"):
        exit("An error occurred: AFLAP_tmp/LA.txt not found. Rerun 01_JELLYFISH.py.")
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
        fa_file = f"AFLAP_tmp/02/{G}_m{args.kmer}_L{LO}_U{UP}.fa"
        if os.path.exists(fa_file):
            print(f"\t\t{args.kmer}-mers for {G} detected. Skipping.")
            histo_same = True
        else:
            print(f"\t\tRunning jellyfish dump for {G}...")
            cmd = f"jellyfish dump -U {UP} -L {LO} -o {fa_file} AFLAP_tmp/01/F0Count/{G}.jf{args.kmer}"
            subprocess.run(cmd, shell=True, executable="/bin/bash")
            histo_same = False

        # counting k-mers
        print(f"\tCounting number of {args.kmer}-mers for {G}...")
        print(f"\t\t{int(os.path.getsize(fa_file) / 2)} {args.kmer}-mers counted for {G}")

        # create histo.png
        png_file = f"AFLAP_Results/Plots/{G}_m{args.kmer}_L{LO}_U{UP}_histo.png"
        if os.path.exists(png_file) and histo_same:
            print(f"\tHistogram for {G} detected. Skipping.")
        else:
            if not histo_same:
                print(f"\tMaking new histogram for {G}...")
            else:
                print(f"\t Making histogram for {G}...")
            histoplot(f"AFLAP_tmp/02/F0Histo/{G}.{args.kmer}.histo", LO, UP, png_file)

            # check if histogram had been built
            if not os.path.exists(png_file) or not os.path.getsize(png_file):
                exit(f"An error occurred: {png_file} not found.")
            print(f"\tHistogram for {G} constructed.")
