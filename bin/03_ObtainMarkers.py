import argparse
import os
import shutil
import subprocess
import sys

def main()->None:
    parser = argparse.ArgumentParser(prog='ObtainMarkers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes.")
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args

    # 1. make directories
    os.makedirs("AFLAP_tmp/03/F0Markers", exist_ok=True)

    # 2. assemble for markers for parents whose bounds are identified
    if not os.path.exists("AFLAP_tmp/01/LA.txt"):
        print("Error in 03_ObtainMarkers.py: AFLAP_tmp/01/LA.txt not found. Rerun 01_JELLYFISH.py.")
    with open("AFLAP_tmp/01/LA.txt", 'r') as fla:
        for p in fla:
            p = p.strip().split()
            G = p[0]    # parent to be evaluated
            LO = p[1]   # lower bound
            UP = p[2]   # upper bound

            # initialize p0 which consists of parents G crosses with
            p0 = []

            # determine what parents G crosses with
            with open("AFLAP_tmp/01/Crosses.txt", 'r') as fc, open(f"AFLAP_tmp/03/{G}_CrossedTo.txt", 'w') as fct:
                for pair in fc:
                    pair = pair.strip().split()

                    if   pair[0] == G: 
                        fct.write(f"{pair[1]}\n")
                        p0.append(pair[1])
                    elif pair[1] == G:
                        fct.write(f"{pair[0]}\n")
                        p0.append(pair[0])
                    else: continue

            # convert p0 into string
            p0 = 'x'.join(p0)

            # check if marker for G exists
            if os.path.exists(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{op}.fa"):
                # check if file is empty
                if not os.path.getsize(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{op}.fa"):
                    print(f"Error in 03_ObtainMarkers.py: AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{op}.fa is empty.")
                    sys.exit(1)

                while True:
                    mark_again = input(f"Marker already found for {G}. Would you like to make a new one? (y/n)")
                    if mark_again in {"yes", "y", "ye"}:
                        print("Performing marker assembly again...")
                        mark_check = True
                    elif mark_again in {"no", "n"}:
                        print("Skipping.")
                        mark_check = False
            else:
                print("Performing marker assembly...")
                mark_check = False

            if mark_check:
                # define ak
                ak = 2 * int(args.kmer) - 1

                # initialize marker stats report variables
                ml_count   = 0      # number of k-mers inputted for assembly
                frag_count = 0      # number of fragments assembled
                frag61     = 0      # number of frags == ak
                frag62     = 0      # number of frags > ak
                mar_count  = 0      # number of markers after refiltering
                mar61      = 0      # number of markers == ak
                mar62      = 0      # number of markers > ak

                # copy .fa file
                if not os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa"):
                    print(f"Error in 03_ObtainMarkers.py: AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa not found. Rerun 02_ExtractSingleCopyMers.py.")
                    sys.exit(1)
                else:
                    shutil.copy(f"AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa", f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa")

                # filter .fa file against other parents
                with open(f"AFLAP_tmp/03.{G}_CrossedTo.txt", 'r') as fct:
                    for op in fct:
                        op = op.strip()

                        # check if .jf for other parent exists
                        if not os.path.exists(f"AFLAP_tmp/01/F0Counts/{op}.jf{args.kmer}"):
                            print(f"Error in 03_ObtainMarkers.py: AFLAP_tmp/01/F0Counts/{op}.jf{args.kmer} not found. Rerun 01_JELLYFISH.py.")
                            sys.exit(1)

                        # filter and overwrite .fa file
                        print(f"\tFiltering and overwriting AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa...")
                        jf_cmd = f"jellyfish query -s AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa AFLAP_tmp/01/F0Counts/{op}.jf{args.kmer}"
                        jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                        with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa", 'w') as ffa:
                            for line in jf_out:
                                line = line.strip().split()

                                # disregard empty lines FIXME: determine why this happens
                                if not len(line): continue

                                # write unique sequences into .fa file
                                if not int(line[0]):
                                    # update ml_count in stats
                                    ffa.write(f">{ml_count}\n{line[0]}\n")
                                    ml_count += 1

                            print(f"\t{ml_count} {G} {args.kmer}-mers remain after filtering against {op}.")

                # initialize k variable
                if   int(args.kmer) == 31: k = 25
                elif int(args.kmer) == 25: k = 19
                else: k = int(args.kmer) - 2

                # perform ABySS assembly
                print(f"\tRunning ABySS with -k set to {k}...")
                abyss_cmd = f"ABYSS -k {k} -c 0 -e 0 AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa -o AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa"
                subprocess.run(abyss_cmd, shell=True)
                # check if abyss ran properly
                if not os.path.exists(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa"):
                    print(f"Error in 03_ObtainMarkers.py: ABySS did not create {G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa.")
                elif not os.path.getsize(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa"):
                    print(f"Error in 03_ObtainMarkers.py: {G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa is empty.")

                # extract subsequences
                print("\tExtracting subsequences...")
                with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa", 'r') as fab, open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss_subseqs.fa", 'w') as fabsub:
                    while True:
                        m = fab.readline()
                        if not m: break

                        if m.startswith('>'):
                            # increment frag stats
                            frag_count += 1

                            # subsequence to abyss subsequence file
                            m = m.strip().replace('>', '').split()
                            if int(m[1]) >= ak:
                                fabsub.write(f">{m[0]}_{m[1]}\n")

                                # define subsequence and its reverse complement
                                subseq = fab.readline().strip()[9:(9 + int(args.kmer))]
                                rc_subseq = subseq[::-1].translate(subseq.maketrans("ATCG", "TAGC"))

                                # compare subsequence and its reverse complement (choose first typographically)
                                if subseq <= rc_subseq: fabsub.write(f"{subseq}\n")
                                else: fabsub.write(f"{rc_subseq}\n")

                                # update frag61 and frag62
                                if int(m[1]) == ak: frag61 += 1
                                else: frag62 += 1

                # refilter against self
                print("\tRefiltering against self...")
                jf_cmd = f"jellyfish query -s AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss_subseqs.fa AFLAP_tmp/01/F0Counts/{G}.jf{args.kmer}"
                jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_jf_query.fa", 'w') as fjq:
                    i = 1
                    for line in jf_out:
                        line = line.strip().split()

                        # disregard empty lines FIXME: determine why this happens
                        if not len(line): continue

                        if int(line[1]) >= int(LO) and int(line[1]) <= int(UP):
                            fjq.write(f">{i}\n{line[0]}\n")
                            i += 1

                # refilter against other parents
                print("\tRefiltering against other parents...")
                with open(f"AFLAP_tmp/03/{G}_CrossedTo.txt", 'r') as fct:
                    for op in fct:
                        op = op.strip()

                        jf_cmd = f"jellyfish query -s AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_jf_query.fa AFLAP_tmp/01/F0Counts/{op}.jf{args.kmer}"
                        jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                        with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_jf_query.fa", 'w') as fjq:
                            i = 1
                            for line in jf_out:
                                line.strip().split()

                                # disregard empty lines FIXME: determine why this happens
                                if not len(line): continue

                                # write unique sequences back into file
                                if not int(line[1]):
                                    fjq.write(f">{i}\n{line[0]}\n")
                                    i += 1
                    
                # create final marker
                print("\tCreating final marker...")
                with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss_subseqs.fa", 'r') as fabsub, open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_jf_query.fa", 'r') as fjq, open(f"AFLAP_tmp/03/ParentalMarkers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{p0}.fa", 'w') as fmark:
                    # create set of sequences from jf_query file
                    fjq_set = set()
                    for line in fjq:
                        if line.startswith('>'): continue
                        fjq_set.add(line.strip())
                    
                    while True:
                        head = fabsub.readline().strip()
                        if not head: break

                        # add sequence to markers file if found in jf_query file
                        seq = fabsub.readline().strip()
                        if seq in fjq_set:
                            fmark.write(f"{head}\n{seq}\n")

                            # update stats
                            mar_count += 1
                            head = head.split('_')

                            if int(head[1]) == ak:
                                mar61 += 1
                            elif int(head[1]) > ak:
                                mar62 += 1

                # print stats
                print(f"Report for {G}:\n" +
                    f"\tNumber of {args.kmer}-mers input into assembly: {ml_count}\n" +
                    f"\tNumber of fragments assembled: {frag_count}\n" +
                    f"\tNumber of fragments == {ak} bp: {frag61}\n" +
                    f"\tNumber of fragments > {ak} bp: {frag62}\n" +
                    f"\tNumber of markers after refiltering: {mar_count}\n" +
                    f"\tNumber of markers == {ak} bp: {mar61}\n" +
                    f"\tNumber of markers > {ak} bp: {mar62}\n")
                
                # write to G's MarkerReport.txt
                with open(f"{G}.MarkerReport.txt", 'w') as f:
                    f.write(f"Report for {G}:\n" +
                            f"\tNumber of {args.kmer}-mers input into assembly: {ml_count}\n" +
                            f"\tNumber of fragments assembled:                  {frag_count}\n" +
                            f"\tNumber of fragments == {ak} bp:                 {frag61}\n" +
                            f"\tNumber of fragments > {ak} bp:                  {frag62}\n" +
                            f"\tNumber of markers after refiltering:            {mar_count}\n" +
                            f"\tNumber of markers == {ak} bp:                   {mar61}\n" +
                            f"\tNumber of markers > {ak} bp:                    {mar62}\n")


if __name__ == "__main__":
    main()