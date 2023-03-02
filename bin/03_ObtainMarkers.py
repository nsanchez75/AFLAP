import argparse
import os
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(prog='ObtainMarkers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes")
    parser.add_argument('-P', '--Pedigree', required=True, help='Pedigree file (required). See AFLAP README for more information.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31]')
    args = parser.parse_args()


    # check for required files from previous scripts
    if not os.path.exists("AFLAP_tmp/01/LA.txt"):
        print("Could not find AFLAP_tmp/01/LA.txt. Please rerun 01_JELLYFISH.py")
        exit(1)

    # dependency check
    try:
        abyss_version = subprocess.check_output("ABYSS --version", shell=True)
        abyss_version = abyss_version.decode().strip()
        print(f"{abyss_version} detected\n")
    except OSError:
        print("ABySS not detected, please modify your PATH. Terminating.")
        exit(1)


    # make new tmp directories
    os.makedirs("AFLAP_tmp/03/ParentalMarkers", exist_ok=True)

    # only assemble markers for those whose boundaries were identified.
    with open("AFLAP_tmp/01/LA.txt", 'r') as fla:
        for g in fla:
            g = g.strip()
            print(f"Beginning analysis for {g}...")

            # initialize stats report variables
            ml_count =   0
            frag_count = 0
            count61 =    0
            count62 =    0
            mar_count =  0
            mar61 =      0
            mar62 =      0

            # identify upper and lower boundaries for file ID and filtering
            print(f"Identifying {args.kmer}-mer boundaries used previously.")
            with open("AFLAP_tmp/02/Boundaries.txt", 'r') as fbound:
                for b in fbound:
                    b = b.strip().split()
                    if b[0] == g:
                        lo = b[1]
                        hi = b[2]
                        break
            print(f"Lower boundary set to {lo}\n" +
                  f"Upper boundary set to {hi}")

            if not os.path.exists(f"AFLAP_tmp/02/ParentalHisto/{g}_m{args.kmer}_L{lo}_U{hi}.fa"):
                print("Cannot find results with these boundaries. Please try to rerun the pipeline. " +
                      "Terminating.")
                exit(1)
            print("Results with these boundaries detected. Proceeding.")

            # identify opposing parental hashes to filter against
            print(f"Identifying parents crossed to {g}...")
            with open("AFLAP_tmp/01/Crosses.txt", 'r') as fin, open(f"AFLAP_tmp/03/{g}_CrossedTo.txt", 'w') as fout:
                for cross in fin:
                    cross = cross.strip().split()
                    if (cross[2] == g and cross[3] == g):
                        print(f"Error: Cannot have {g} as both parents in AFLAP_tmp/01/Crosses.txt." +
                               "Terminating.")
                        exit(1)
                    elif (cross[2] == g): fout.write(f"{cross[3]}\n")
                    elif (cross[3] == g): fout.write(f"{cross[2]}\n")
                    else: continue

            # copy .fa file into 03 for rewriting purposes
            shutil.copy(f"AFLAP_tmp/02/ParentalHisto/{g}_m{args.kmer}_L{lo}_U{hi}.fa", f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}.fa")

            # filter against opposing parents
            with open(f"AFLAP_tmp/03/{g}_CrossedTo.txt", 'r') as f:
                for p0 in f:
                    p0 = p0.strip()

                    if os.path.exists(f"AFLAP_tmp/03/ParentalMarkers/{g}_m{args.kmer}_MARKERS_L{lo}_U{hi}_{p0}.fa"):
                        print("Previously calculated markers detected. If you want to calculate new markers, please delete:\t" +
                             f"./AFLAP_tmp/03/ParentalMarkers/{g}_m{args.kmer}_MARKERS_L{lo}_U{hi}_{p0}.fa\n" +
                             f"Summary of markers available:\n\t./{g}.MarkerReport.txt")
                    else:
                        if not os.path.exists(f"AFLAP_tmp/01/ParentalCounts/{p0}.jf{args.kmer}"):
                            print(f"JELLYFISH results of {p0} can't be located. Please rerun 01_JELLYFISH.py. " +
                                   "Terminating.")
                            exit(0)
                        print(f"Intersecting {g} with {p0}...")

                        # filter and overwrite output
                        print("Filtering and overwriting output...")
                        jellyfish_cmd = f"jellyfish query -s AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}.fa AFLAP_tmp/01/ParentalCounts/{p0}.jf{args.kmer}"
                        jellyfish_output = subprocess.run(jellyfish_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                        with open(f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}.fa", 'w') as fout:
                            for line in jellyfish_output:
                                line = line.strip().split()

                                # disregard empty list at end of jellyfish_output FIXME: figure out why this is happening
                                if len(line) == 0:
                                    continue

                                if int(line[1]) == 0:
                                    fout.write(f">{ml_count}\n{line[0]}\n")
                                    ml_count += 1
                            print(f"{ml_count} {g} {args.kmer}-mers remain after filtering against {p0}.")

            # assemble using ABySS (optimal k for 31-mer is 25 and For 21-mer is 19) // mark1
            if int(args.kmer) == 31:   k = 25
            elif int(args.kmer) == 25: k = 19
            else:                 k = int(args.kmer) - 2

            print(f"Running ABySS with -k set to {k}...")
            abyss_cmd = f"ABYSS -k {k} -c 0 -e 0 AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}.fa -o AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark1.fa"
            subprocess.run(abyss_cmd)
            abyss_check = False

            # extract subsequences // mark2
            ak = int(args.kmer) + int(args.kmer) - 1
            print("Extracting subsequences...")
            with open(f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark1.fa", 'r') as fmark1, open(f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark2.fa", 'w') as fmark2:
                m = fmark1.readline()
                while m:
                    if m.startswith('>'):
                        m = m.strip().replace('>', '').split()
                        if int(m[1]) >= ak:
                            fmark2.write(f">{m[0]}_{m[1]}\n")
                            seq = fmark1.readline()
                            subseq = seq[9:(9 + int(args.kmer))]

                            # find the reverse complement and write to file whichever one is lexicographically first
                            rc_subseq = subseq[::-1].translate(subseq.maketrans("ATCG", "TAGC"))

                            if subseq <= rc_subseq:
                                fmark2.write(f"{subseq}\n")
                            else:
                                fmark2.write(f"{rc_subseq}\n")

                        # modify stats
                        frag_count += 1
                        if int(m[1]) == ak:
                            count61 += 1
                        elif int(m[1]) > ak:
                            count62 += 1

                    m = fmark1.readline()

            # refilter against self for k-mers between boundaries using jellyfish // mark3
            print("Refiltering against self...")
            jellyfish_cmd = f"jellyfish query -s AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark2.fa AFLAP_tmp/01/ParentalCounts/{g}.jf{args.kmer}"  # FIXME: only works if jellyfish is topological
            jellyfish_output = subprocess.run(jellyfish_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
            with open(f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark3.fa", 'w') as fmark3:
                i = 1
                for line in jellyfish_output:
                    line = line.strip().split()

                    # disregard empty list at end of jellyfish_output FIXME: figure out why this is happening
                    if len(line) == 0:
                        continue

                    if int(line[1]) >= int(lo) and int(line[1]) <= int(hi):

                        fmark3.write(f">{i}\n{line[0]}\n")
                        i += 1

            # refilter against other parents using jellyfish // mark3
            print("Refiltering againts other parents...")
            with open(f"AFLAP_tmp/03/{g}_CrossedTo.txt", 'r') as fin:
                for p in fin:
                    p = p.strip()

                    jellyfish_cmd = f"jellyfish query -s AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark3.fa AFLAP_tmp/01/ParentalCounts/{p}.jf{args.kmer}"
                    jellyfish_output = subprocess.run(jellyfish_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')

                    with open(f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark3.fa", 'w') as fmark3:
                        i = 1
                        for line in jellyfish_output:
                            line = line.strip().split()

                            if len(line) >= 2 and int(line[1]) == 0:
                                fmark3.write(f">{i}\n{line[0]}\n")

            # export final marker set with ABySS-conserved headers
            with open(f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark2.fa", 'r') as fmark2, open(f"AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark3.fa", 'r') as fmark3, open(f"AFLAP_tmp/03/ParentalMarkers/{g}_m{args.kmer}_MARKERS_L{lo}_U{hi}_{p0}.fa", 'w') as fout:
                # create set of Mark3 file sequences
                fmark3_set = set()
                for line in fmark3:
                    if line.startswith('>'):
                        continue
                    line = line.strip()
                    fmark3_set.add(line)

                head = fmark2.readline().strip()
                while head:
                    seq = fmark2.readline().strip()

                    # add sequence to MARKERS file if found in both Mark2 and Mark3 FIXME: MARKERS doesn't have anything
                    if seq in fmark3_set:
                        fout.write(f"{head}\n{seq}\n")

                        # modify stats
                        mar_count += 1
                        head = head.split('_')

                        if int(head[1]) == ak:
                            mar61 += 1
                        elif int(head[1]) > ak:
                            mar62 += 1

                    head = fmark2.readline().strip()

            # print stats
            if abyss_check:
                print(f"Report for {g}:\n" +
                    f"\tNumber of {args.kmer}-mers input into assembly: {ml_count}\n" +
                    f"\tNumber of fragments assembled: {frag_count}\n" +
                    f"\tNumber of fragments == {ak} bp: {count61}\n" +
                    f"\tNumber of fragments > {ak} bp: {count62}\n" +
                    f"\tNumber of markers after refiltering: {mar_count}\n" +
                    f"\tNumber of markers == {ak} bp: {mar61}\n" +
                    f"\tNumber of markers > {ak} bp: {mar62}\n")

                # write to MarkerReport.txt
                with open(f"{g}.MarkerReport.txt", 'w') as f:
                    f.write(f"Report for {g}:\n" +
                            f"\tNumber of {args.kmer}-mers input into assembly: {ml_count}\n" +
                            f"\tNumber of fragments assembled:                  {frag_count}\n" +
                            f"\tNumber of fragments == {ak} bp:                 {count61}\n" +
                            f"\tNumber of fragments > {ak} bp:                  {count62}\n" +
                            f"\tNumber of markers after refiltering:            {mar_count}\n" +
                            f"\tNumber of markers == {ak} bp:                   {mar61}\n" +
                            f"\tNumber of markers > {ak} bp:                    {mar62}\n")
            else:
                print(f"ABySS of {g} was detected. " +
                      f"If you would like to create a new marker report, delete AFLAP_tmp/03/{g}_m{args.kmer}_L{lo}_U{hi}_Mark1.fa.")

if __name__ == "__main__":
    main()