import argparse
import multiprocessing as mp
import os
import pandas as pd
import shutil
import subprocess

from get_LA_info import get_LA_info
from find_identical_loci import find_identical_loci

#################################################
#	A Python script to derive single copy k-mers that are unique to a parent. These are then used a markers.
#	To reduce redundancy k-mers are assembled using ABySS.
#	To enable the use of a consistent has size, the markers are reduced to a sequnce length equal to option m.
#	This consistent marker length then only need to be surveyed against only one progeny hash.
#################################################

def abyss_assembly(k:int, G:str, LO:int, UP:int, kmer:int, abyss_file:str, fafile:str)->None:
    # check if abyss file for G already exists
    if os.path.exists(abyss_file) and os.path.getsize(abyss_file):
        print(f"ABySS files for {G} detected. Skipping...")
    # run abyss otherwise
    else:
        print(f"\tRunning ABySS with -k set to {k}...\n")
        subprocess.run(args=f"ABYSS -k {k} -c 0 -e 0 {fafile} -o AFLAP_tmp/03/{G}_m{kmer}_L{LO}_U{UP}_abyss.fa",
                    shell=True)
        # check if abyss ran properly
        if not os.path.exists(abyss_file):
            exit(f"An error occurred: ABySS did not create {abyss_file}.")
        elif not os.path.getsize(abyss_file):
            exit(f"An error occurred: {abyss_file} is empty.")

def get_markers(G_info:tuple, kmer:int)->None:
    G, LO, UP, P0 = G_info
    seq_groups = pd.DataFrame(columns=["Sequence", "Locus Sequence"])
    ak = 2 * int(kmer) - 1

    print(f"Performing marker assembly on {G}...")

    # initialize marker stats report variables
    ml_count = frag_count = frag61 = frag62 = mar_count = mar61 = mar62 = 0

    # copy .fa file
    fafile_02 = f"AFLAP_tmp/02/{G}_m{kmer}_L{LO}_U{UP}.fa"
    fafile_03 = f"AFLAP_tmp/03/{G}_m{kmer}_L{LO}_U{UP}.fa"
    if not os.path.exists(fafile_02):
        exit(f"An error occurred: {fafile_02} not found. Rerun 02_ExtractSingleCopyMers.py.")
    shutil.copy(fafile_02, fafile_03)

    # filter .fa file against other parents
    for op in P0.split():
        op = op.strip()

        # check if .jf for other parent exists
        opfile = f"AFLAP_tmp/01/F0Count/{op}.jf{kmer}"
        if not os.path.exists(opfile):
            exit(f"An error occurred: {opfile} not found. Rerun 01_JELLYFISH.py.")

        # filter and overwrite .fa file
        print(f"\tFiltering and overwriting {fafile_03}...")
        jf_out = subprocess.run(args=f"jellyfish query -s {fafile_03} {opfile}",
                                shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
        with open(f"{fafile_03}", 'w') as ffa:
            for line in jf_out:
                line = line.strip().split()

                # write unique sequences into .fa file
                if len(line) and not int(line[1]):
                    # update ml_count in stats
                    ffa.write(f">{ml_count}\n{line[0]}\n")
                    ml_count += 1

            print(f"\t{ml_count} {G} {kmer}-mers remain after filtering against {op}.")

    # perform ABySS assembly
    if   int(kmer) == 31: k = 25
    elif int(kmer) == 25: k = 19
    else: k = int(kmer) - 2
    abyss_file = f"AFLAP_tmp/03/{G}_m{kmer}_L{LO}_U{UP}_abyss.fa"
    abyss_assembly(k, G, LO, UP, kmer, abyss_file, fafile_03)

    # extract subsequences
    print("\tExtracting subsequences...")
    abyss_subseq_file = f"AFLAP_tmp/03/{G}_m{kmer}_L{LO}_U{UP}_abyss_subseqs.fa"
    with open(abyss_file, 'r') as fab, open(abyss_subseq_file, 'w') as fabsub:
        while True:
            id = fab.readline().strip()
            seq = fab.readline().strip()
            if not id or not seq: break

            # increment frag stats
            frag_count += 1

            # analyze sequence if passes against ak
            id = id.replace('>', '').split()
            if int(id[1]) >= ak:
                # get id
                fabsub.write(f">{id[0]}_{id[1]}\n")

                # get subsequence and its reverse complement
                subseq = seq[9:(9 + int(kmer))]
                rc_subseq = subseq[::-1].translate(subseq.maketrans("ATCG", "TAGC"))
                # choose subsequence by first alphabetically
                if subseq > rc_subseq: subseq = rc_subseq
                fabsub.write(f"{subseq}\n")

                # get sequence locus via first and last couple of base pairs
                seq_groups.loc[len(seq_groups.index)] = [subseq, seq[0:(int(kmer) - 1)] + seq[(len(seq) - int(kmer) + 1):]]

                # update frag61 and frag62
                if int(id[1]) == ak: frag61 += 1
                else: frag62 += 1

    # refilter against self
    print("\tRefiltering against self...")
    jf_cmd = f"jellyfish query -s {abyss_subseq_file} AFLAP_tmp/01/F0Count/{G}.jf{kmer}"
    jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
    jfqfile = f"AFLAP_tmp/03/{G}_m{kmer}_L{LO}_U{UP}_jf_query.fa"
    with open(jfqfile, 'w') as fjq:
        i = 1
        for line in jf_out:
            line = line.strip().split()

            if len(line) and int(line[1]) >= int(LO) and int(line[1]) <= int(UP):
                fjq.write(f">{i}\n{line[0]}\n")
                i += 1

    # refilter against other parents
    print("\tRefiltering against other parents...")
    for op in P0.split():
        op = op.strip()

        jf_cmd = f"jellyfish query -s {jfqfile} AFLAP_tmp/01/F0Count/{op}.jf{kmer}"
        jf_out = subprocess.run(jf_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
        with open(jfqfile, 'w') as fjq:
            i = 1
            for line in jf_out:
                line = line.strip().split()

                # write unique sequences back into file
                if len(line) and not int(line[1]):
                    fjq.write(f">{i}\n{line[0]}\n")
                    i += 1

    # create final marker
    print("\tCreating final marker...")
    with open(abyss_subseq_file, 'r') as fabsub, open(jfqfile, 'r') as fjq, open(f"AFLAP_tmp/03/F0Markers/{G}_m{kmer}_MARKERS_L{LO}_U{UP}_{P0}.fa", 'w') as fmark:
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

    stats = f"Report for {G}:\n" + \
            f"\tNumber of {kmer}-mers input into assembly:      {ml_count}\n" + \
            f"\tNumber of fragments assembled:                  {frag_count}\n" + \
            f"\tNumber of fragments == {ak} bp:                 {frag61}\n" + \
            f"\tNumber of fragments > {ak} bp:                  {frag62}\n" + \
            f"\tNumber of markers after refiltering:            {mar_count}\n" + \
            f"\tNumber of markers == {ak} bp:                   {mar61}\n" + \
            f"\tNumber of markers > {ak} bp:                    {mar62}\n"

    # print and write stats to parent's MarkerReport.txt
    print(stats)
    with open(f"{G}.MarkerReport.txt", 'w') as f: f.write(stats)

    # convert seq_groups dataframe to txt file
    seq_groups.to_csv(f"AFLAP_tmp/04/SimGroups/{G}_locus_seqs.txt", sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ObtainMarkers', description="A script to obtain single copy k-mers from parental JELLYFISH hashes.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()

    # make directories
    os.makedirs("AFLAP_tmp/03/F0Markers", exist_ok=True)
    os.makedirs("AFLAP_tmp/03/SimGroups", exist_ok=True)

    # assemble for markers for parents whose bounds are identified
    processes = list()
    for G_info in get_LA_info():
        p = mp.Process(target=get_markers, args=(G_info, args.kmer))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
    print("Identified markers for all parents.")

    # TODO: find sequences of identical loci