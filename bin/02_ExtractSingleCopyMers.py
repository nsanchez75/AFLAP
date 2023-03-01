import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(prog='ExtractSingleCopyMers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes.")
    parser.add_argument('-P', '--Pedigree', required=True, help='Pedigree file (required). See AFLAP README for more information.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31]')
    args = parser.parse_args()


    # 1. Make new directories for the second script
    os.makedirs("AFLAP_Results/Plots", exist_ok=True)
    os.makedirs("AFLAP_tmp/02/ParentalHisto", exist_ok=True)


    # 2. Obtain histograms
    print("Generating histograms for F0 to undergo linkage analysis...")

    # with open("AFLAP_tmp/01/Crosses.txt", 'r') as f:
    #     lines = f.readlines()

    #     f1_values = set()
    #     f2_values = set()
    #     for line in lines:
    #         line_split = line.split()
    #         if line_split[1] == '1':
    #             f1_values.add((line_split[2], line_split[3]))
    #         elif line_split[1] == '2':
    #             f2_values.add((line_split[2], line_split[3]))

    #     f1 = len(f1_values)
    #     f2 = len(f2_values)

    with open("AFLAP_tmp/01/LA.txt") as f:
        for line in f:
            g = line.strip()
            # check if it has already been created
            if os.path.exists(f"AFLAP_tmp/02/ParentalHisto/{g}.{args.kmer}.histo"):     # TODO: fix skip check
                print(f"Histogram for {g} {args.kmer}-mer detected. Skipping." +
                      f"\n\tNot correct? Cancel and delete AFLAP_tmp/02/ParentalHisto/{g}.{args.kmer}.histo or rerun in a clean directory.")
            else:
                cmd = f"jellyfish histo AFLAP_tmp/01/ParentalCounts/{g}.jf{args.kmer} > AFLAP_tmp/02/ParentalHisto/{g}.{args.kmer}.histo"
                subprocess.run(cmd, shell=True, executable="/bin/bash")
                print(f"Histogram for {g} generated.")

            # check to see if coverage boundaries provided by user
            with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as f:
                for parent in f:
                    parent = parent.strip().split()
                    if len(parent) == 5 and parent[0] == g:
                        print(f"User has supplied lower and upper boundaries for {g} {args.kmer}-mer cut-offs. On to extraction.")
                        lo = int(parent[3])
                        hi = int(parent[4])
                    else:
                        # where super hacky peak finding and boundary setting would've been
                        continue

            print(f"\nLower boundary for {g} set to {lo}, upper boundary to {hi}.")
            with open("AFLAP_tmp/02/Boundaries.txt", 'w') as f:
                f.write(f"{g}\t{lo}\t{hi}")

            print(f"Extracting {args.kmer}-mers from {g}...")

            if os.path.exists(f"AFLAP_tmp/02/ParentalHisto/{g}_m{args.kmer}_L{lo}_U{hi}.fa"):
                print(f"{g} {args.kmer}-mer previously extracted between {lo} and {hi}.\n" +
                      f"\tDelete {g}_m{args.kmer}_L{lo}_U{hi}.fa to rebuild.")
            else:
                cmd = f"jellyfish dump -U {hi} -L {lo} -o AFLAP_tmp/02/ParentalHisto/{g}_m{args.kmer}_L{lo}_U{hi}.fa AFLAP_tmp/01/ParentalCounts/{g}.jf{args.kmer}"
                subprocess.run(cmd, shell=True, executable="/bin/bash")

            mer_count = 0
            with open(f"AFLAP_tmp/02/ParentalHisto/{g}_m{args.kmer}_L{lo}_U{hi}.fa", 'r') as f:
                for mer in f:
                    if mer.startswith('>'): mer_count += 1
            print(f"\t{mer_count} {args.kmer}-mers extracted from {g}.\n")

            if os.path.exists(f"AFLAP_Results/Plots/{g}_m{args.kmer}_L{lo}_U{hi}_histo.png"):
                make_histo = input(f"A histogram has been detected for {g} with Lower={lo} and Upper={hi}. Would you like to make a new one? (y/n) ")
                while True:
                    if make_histo.lower() in ['n', 'no']:
                        break
                    elif make_histo.lower() in ['y', 'yes']:
                        print("Producing parental histogram...")
                        cmd = f"Rscript bin/HistoPlot.R AFLAP_tmp/02/ParentalHisto/{g}.{args.kmer}.histo {lo} {hi} AFLAP_Results/Plots/{g}_m{args.kmer}_L{lo}_U{hi}_histo.png"
                        subprocess.run(cmd, shell=True)
                        break
                    else:
                        print("Invalid input. Type (y/n) ")
            else:
                print("Producing parental histogram...")
                cmd = f"Rscript bin/HistoPlot.R AFLAP_tmp/02/ParentalHisto/{g}.{args.kmer}.histo {lo} {hi} AFLAP_Results/Plots/{g}_m{args.kmer}_L{lo}_U{hi}_histo.png"
                subprocess.run(cmd, shell=True)


if __name__ == "__main__":
    main()