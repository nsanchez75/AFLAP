import argparse
import os
import pandas as pd
import shutil
import subprocess

from get_LA_info import get_LA_info

#################################################
#	A Python script to derive single copy k-mers that are unique to a parent. These are then used a markers.
#	To reduce redundancy k-mers are assembled using ABySS.
#	To enable the use of a consistent has size, the markers are reduced to a sequnce length equal to option m.
#	This consistent marker length then only need to be surveyed against only one progeny hash.
#################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ObtainMarkers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()

    # make directories
    os.makedirs("AFLAP_tmp/03/F0Markers", exist_ok=True)
    os.makedirs("AFLAP_tmp/03/SimGroups", exist_ok=True)

    # assemble for markers for parents whose bounds are identified
    try:
        # initialize a sequence grouper dataframe
        seq_groups = pd.DataFrame(columns=["Sequence", "Parent", "Identifier"])

        for G_info in get_LA_info():
            G, LO, UP, P0 = G_info

            print(f"Performing marker assembly on {G}...")
            # define ak
            ak = 2 * int(args.kmer) - 1

            # initialize marker stats report variables
            ml_count   = 0  # number of k-mers inputted for assembly
            frag_count = 0  # number of fragments assembled
            frag61     = 0  # number of frags == ak
            frag62     = 0  # number of frags > ak
            mar_count  = 0  # number of markers after refiltering
            mar61      = 0  # number of markers == ak
            mar62      = 0  # number of markers > ak

            # copy .fa file
            if not os.path.exists(f"AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa"):
                raise FileNotFoundError(f"AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa not found. Rerun 02_ExtractSingleCopyMers.py.")
            shutil.copy(f"AFLAP_tmp/02/F0Histo/{G}_m{args.kmer}_L{LO}_U{UP}.fa",
                        f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa")

            # filter .fa file against other parents
            for op in P0.split():
                op = op.strip()

                # check if .jf for other parent exists
                if not os.path.exists(f"AFLAP_tmp/01/F0Count/{op}.jf{args.kmer}"):
                    raise FileNotFoundError(f"AFLAP_tmp/01/F0Count/{op}.jf{args.kmer} not found. Rerun 01_JELLYFISH.py.")

                # filter and overwrite .fa file
                print(f"\tFiltering and overwriting AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa...")
                jf_cmd = f"jellyfish query -s AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa AFLAP_tmp/01/F0Count/{op}.jf{args.kmer}"
                jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa", 'w') as ffa:
                    for line in jf_out:
                        line = line.strip().split()

                        # disregard empty lines FIXME: determine why this happens
                        if not len(line): continue

                        # write unique sequences into .fa file
                        if not int(line[1]):
                            # update ml_count in stats
                            ffa.write(f">{ml_count}\n{line[0]}\n")
                            ml_count += 1

                    print(f"\t{ml_count} {G} {args.kmer}-mers remain after filtering against {op}.")

            # initialize k variable
            if   int(args.kmer) == 31: k = 25
            elif int(args.kmer) == 25: k = 19
            else: k = int(args.kmer) - 2

            # perform ABySS assembly
            print(f"\tRunning ABySS with -k set to {k}...\n")
            abyss_cmd = f"ABYSS -k {k} -c 0 -e 0 AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}.fa -o AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa"
            subprocess.run(abyss_cmd, shell=True)
            # check if abyss ran properly
            if not os.path.exists(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa"):
                raise FileNotFoundError(f"ABySS did not create {G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa.")
            elif not os.path.getsize(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa"):
                raise ValueError(f"{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa is empty.")

            # extract subsequences
            print("\tExtracting subsequences...")
            with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss.fa", 'r') as fab, open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss_subseqs.fa", 'w') as fabsub:
                while True:
                    m = fab.readline()
                    if not m: break
                    if not m.startswith('>'): continue

                    # increment frag stats
                    frag_count += 1

                    # get sequence
                    seq = fab.readline().strip()

                    # add (k-mer - 1) of each end of sequence to sequence group dataframe
                    seq_groups.loc[len(seq_groups.index)] = [seq[0:(int(args.kmer) - 1)] + seq[(len(seq) - int(args.kmer) + 1):], G, f">{m[0]}_{m[1]}"]

                    # subsequence to abyss subsequence file
                    m = m.strip().replace('>', '').split()
                    if int(m[1]) >= ak:
                        fabsub.write(f">{m[0]}_{m[1]}\n")

                        # define subsequence and its reverse complement
                        subseq = seq[9:(9 + int(args.kmer))]
                        rc_subseq = subseq[::-1].translate(subseq.maketrans("ATCG", "TAGC"))

                        # compare subsequence and its reverse complement (choose first typographically)
                        if subseq <= rc_subseq: fabsub.write(f"{subseq}\n")
                        else: fabsub.write(f"{rc_subseq}\n")

                        # update frag61 and frag62
                        if int(m[1]) == ak: frag61 += 1
                        else: frag62 += 1

            # refilter against self
            print("\tRefiltering against self...")
            jf_cmd = f"jellyfish query -s AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss_subseqs.fa AFLAP_tmp/01/F0Count/{G}.jf{args.kmer}"
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
            for op in P0.split():
                op = op.strip()

                jf_cmd = f"jellyfish query -s AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_jf_query.fa AFLAP_tmp/01/F0Count/{op}.jf{args.kmer}"
                jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_jf_query.fa", 'w') as fjq:
                    i = 1
                    for line in jf_out:
                        line = line.strip().split()

                        # disregard empty lines FIXME: determine why this happens
                        if not len(line): continue

                        # write unique sequences back into file
                        if not int(line[1]):
                            fjq.write(f">{i}\n{line[0]}\n")
                            i += 1

            # create final marker
            print("\tCreating final marker...")
            with open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_abyss_subseqs.fa", 'r') as fabsub, open(f"AFLAP_tmp/03/{G}_m{args.kmer}_L{LO}_U{UP}_jf_query.fa", 'r') as fjq, open(f"AFLAP_tmp/03/F0Markers/{G}_m{args.kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa", 'w') as fmark:
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

        # # find identical loci
        # print("Finding homozygous sequences...")
        # ## get parents categorized by sex
        # parents_lists = list()
        # with open("AFLAP_tmp/Crosses.txt", 'r') as fcrosses:
        #     for cross in fcrosses:
        #         cross = cross.strip().split()
        #         parents_lists.append([cross[2], cross[3]])
        # ## create file of sequences
        # unique_seqs = pd.DataFrame(columns=["Male Identifier", "Female Identifier", "Sequence"])
        # seq_groups = seq_groups[seq_groups.duplicated("Sequence", keep=False)]
        # for useq in seq_groups["Sequence"].unique():
        #     useq_df = seq_groups[seq_groups["Sequence"] == useq]
        #     for plist in parents_lists:
        #         if set(useq_df["Parent"].unique()) == set(plist):
        #             male_identifier = useq_df[useq_df["Parent"] == plist[0]]["Identifier"]
        #             female_identifier = useq_df[useq_df["Parent"] == plist[1]]["Identifier"]
        #             unique_seqs.loc[len(unique_seqs.index)] = [, useq]
        # unique_seqs.to_csv("AFLAP_tmp/03/Locus.tsv", sep='\t', index=False)
        # print(f"{len(unique_seqs.index)} common loci found.")

    except Exception as e:
        print(f"Error in 03_ObtainMarkers.py: {e}")
        exit(1)