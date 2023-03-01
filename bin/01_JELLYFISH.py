import argparse
import gzip
import os
import subprocess


###########################################################
#	Shell script to be ran as part of the AFLAP pipeline.
#	This script will sanity check the pedigree file then run JELLYFISH on parents and progeny specified in the pedigree file.
#	The output will be stored in AFLAP_tmp.
#	The script will detect previous results and used them when able.
###########################################################

def main():
    parser = argparse.ArgumentParser(prog='JELLYFISH', description='A script to obtain JELLYFISH hashes to be used when running AFLAP.')
    parser.add_argument('-P', '--Pedigree', required=True, help='Pedigree file (required). See AFLAP README for more information.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31]')
    parser.add_argument('-t', '--threads', default=4, help='Threads for JELLYFISH counting (optional). Default [4]')
    args = parser.parse_args()


    # 1. Identify which parents the map is to be constructed of (this info will be used in later scripts for the most part)
    with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as fin, open("AFLAP_tmp/01/NoLA.txt", 'w') as fnola, open("AFLAP_tmp/01/LA.txt", 'w') as fla:
        parents = set()
        for line in fin:
            col = line.strip().split()
            if col[0] not in parents:
                if col[3] == "NA" or col[4] == "NA":
                    fnola.write(f"{col[0]}\n")
                else:
                    fla.write(f"{col[0]}\n")

                parents.add(col[0])

    with open("AFLAP_tmp/01/LA.txt", 'r') as f1, open("AFLAP_tmp/01/NoLA.txt", 'r') as f2:
        lmb = len(f1.readlines())
        nlm = len(f2.readlines())

        if not lmb:
            print("Hmm it seems the user has requested no linkage maps be built with parental derived markers.\n" +
                  "I reccomend you check the pedgree file, at least one F0 should have no NA in fields 4 and 5.\n" +
                  "They can be blank for automatic coverage calculation, or contain coverage cut-off boundaries.\n" +
                  "Terminating.")
            exit(1)

        if nlm:
            print("Linkage analysis will not be performed on these F0:\n")
            f2.seek(0)
            for line in f2:
                line = line.strip()
                print(f"\t\t{line}")
            print(f"\nDoesn't sound right? Make sure 'NA' does not occur in fields 4 and 5 for these F0 in {args.Pedigree}")
        else:
            print(f"We are in business! Linkage analysis will be performed on markers derived from {lmb} F0:")
            f1.seek(0)
            for line in f1:
                print(line)


    # 2. K-Mer counting
    print("\nMoving on to k-mer counting.\n")

    ## make directories
    os.makedirs("AFLAP_tmp/01/ParentalCounts", exist_ok=True)
    os.makedirs("AFLAP_tmp/01/ProgCounts", exist_ok=True)

    ## Parents
    with open("AFLAP_tmp/01/F0.txt", 'r') as f:
        lines = f.readlines()
        pla = len(lines)
        print(f"Beginning k-mer counting for {pla} parents...")

        for line in lines:
            g = line.strip()
            # check if it has already been created
            if os.path.exists(f"AFLAP_tmp/01/ParentalCounts/{g}.jf{args.kmer}"):
                # if exists then skip
                print(f"{args.kmer} hash detected for {g}. Skipping.\n" +
                      f"\tNot correct? Cancel and delete AFLAP_tmp/01/ParentalCounts/{g}.jf{args.kmer} or run from clean directory.")
            # if it doesn't exist then generate
            else:
                print(f"Beginning k-mer counting for {g}.")
                reads = []
                with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as fin:
                    for line in fin:
                        col = line.strip().split()
                        if col[0] == g: reads.append(col[2])

                for r in reads:
                    # check if file exists
                    if not os.path.exists(r):
                        print(f"File {r} not found. Terminating.")
                        exit(1)

                # convert reads to string with spaces between file paths
                j_fin = ' '.join(reads)

                # run jellyfish cound command
                cmd = f"jellyfish count -m {args.kmer} -C -s 1G -t {args.threads} -o AFLAP_tmp/01/ParentalCounts/{g}.jf{args.kmer} <(zcat {j_fin})"
                # print(cmd)      # used to debug jellyfish
                subprocess.run(cmd, shell=True, executable="/bin/bash")

                if not os.path.exists(f"AFLAP_tmp/01/ParentalCounts/{g}.jf{args.kmer}"):
                    print(f"JELLYFISH for {g} did not complete. Is the file gzipped? Terminating.")
                    exit(1)


    ## Progeny
    print("\nParental K-mer counting complete!\nOn to the progeny!\n")

    ### F1 generation
    #### generate F1 progeny file
    with open("AFLAP_tmp/Pedigree_F1.txt") as fin, open("AFLAP_tmp/01/Prog_F1.txt", 'w') as fout:
        prog_f1 = []
        for line in fin:
            col = line.strip().split()
            prog_f1 += [col[0]]

        prog_f1.sort() # FIXME: remove redundant lines
        prog_set = set()
        for i in prog_f1:
            if i not in prog_set:
                fout.write(f"{i}\n")
                prog_set.add(i)

    #### perform k-mer counting on F1 progeny
    with open("AFLAP_tmp/01/Prog_F1.txt", 'r') as f:
        lines = f.readlines()
        prog_count = len(lines)
        print(f"Beginning k-mer counting for {prog_count} F1 progeny...")

        for line in lines:
            g = line.strip()
            # check if it has already been created
            if os.path.exists(f"AFLAP_tmp/01/ProgCounts/{g}.jf{args.kmer}"):
                # skip if true
                print(f"{args.kmer} hash detected for {g}. Skipping.\n" +
                      f"\tNot correct? Cancel and delete AFLAP_tmp/01/ProgCounts/{g}.jf{args.kmer}, or run from clean directory.")
            # generate if false
            else:
                print(f"Beginning k-mer counting for {g}.")
                reads = []
                with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as fin:
                    for line in fin:
                        col = line.strip().split()
                        if col[0] == g: reads.append(col[2])

                for r in reads:
                    if not os.path.exists(r):
                        print(f"File {r} not found. Terminating.")
                        exit(1)

                # convert reads to string with spaces between file paths
                j_fin = ' '.join(reads)

                # run jellyfish count command
                cmd = f"jellyfish count -m {args.kmer} -C -s 1G -t {args.threads} -o AFLAP_tmp/01/ProgCounts/{g}.jf{args.kmer} <(zcat {j_fin})"
                # print(cmd)    # used to debug jellyfish
                subprocess.run(cmd, shell=True, executable="/bin/bash")

                if not os.path.exists(f"AFLAP_tmp/01/ProgCounts/{g}.jf{args.kmer}"):
                    print(f"JELLYFISH for {g} did not complete. Is the file gzipped? Terminating.")
                    exit(1)

    ### F2 generation (TODO)

    print("\nk-mer counting done!")

if __name__ == "__main__":
    main()
