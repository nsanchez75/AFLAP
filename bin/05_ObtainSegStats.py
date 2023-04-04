import argparse
import numpy as np
import os
import re
import subprocess
import sys

def histo_sort(line:str)->float:
    line_fields = line.strip().split()
    return float(line_fields[0])

def main()->None:
    parser = argparse.ArgumentParser(prog='ObtainSegStats', description='A script to plot marker distributions in progeny.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()


    # 1. create directories
    os.makedirs("AFLAP_tmp/05/FilteredCall", exist_ok=True)
    os.makedirs("AFLAP_tmp/05/SegregationInfo", exist_ok=True)


    # 2. initialize variables
    ak = 2 * int(args.kmer) - 1


    # 3. analyze all parents whom we can create a genetic map for
    print("Performing segment statistics analysis...")
    with open("AFLAP_tmp/01/LA.txt", 'r') as fla:
        for p in fla:
            # extract info about analyzed parent
            p = p.strip().split()
            G = p[0]
            LO = p[1]
            UP = p[2]

            print(f"\tObtaining segment statistics for {G}...")

            with open(f"AFLAP_tmp/03/{G}_CrossedTo.txt", 'r') as fct:
                p0 = []
                for op in fct:
                    p0.append(op.strip())
                p0 = '_'.join(p0)

            # check for genotype table
            if not os.path.exists(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.MarkerID.tsv"):
                raise FileNotFoundError(f"Error in 05_ObtainSegStats.py: Genotype table for {G} not found. Rerun 04_Genotyping.py.")
                sys.exit(1)

            # count progeny
            with open(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.MarkerID.tsv", 'r') as ftsv:
                num_prog = len(ftsv.readline()) - 2

            print(f"\t\t{num_prog} Genotype calls for {G} detected. Summarizing...")

            # compare contig lengths that subsequences came from TODO: fact check whether it comes from ABySS or jellyfish
            with open(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.MarkerID.tsv", 'r') as ftsv:
                # initialize dictionaries for MarkerEqual and MarkerOver
                meq = dict()
                mov = dict()


                for tsv_line in ftsv:
                    # extract original contig length
                    tsv_line = tsv_line.strip()
                    # replace underscore with space
                    tsv_line = re.sub('_', ' ', tsv_line)
                    # split lines by whitespace
                    tsv_line = tsv_line.split()

                    counts = sum([int(x) for x in tsv_line[4:]])
                    prop = counts / num_prog

                    # write to MarkerEqual and MarkerOver
                    if int(tsv_line[2]) == ak:
                        if prop in meq:
                            meq[prop] = counts
                        
                        with open(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerEqual{ak}.hist", 'a') as fme:
                            fme.write(f"{prop} {counts}\n")
                    else:
                        with open(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerOver{ak}.hist", 'a') as fmo:
                            fmo.write(f"{prop} {counts}\n")

            # sort histograms (helps with debugging)
            with open(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerEqual{ak}.hist", 'r+') as fme:
                lines = fme.readlines()
                lines.sort(key=histo_sort)
                fme.seek(0)
                fme.writelines(lines)
            with open(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerOver{ak}.hist", 'r+') as fmo:
                fmo.seek(0)
                lines = fmo.readlines()
                lines.sort(key=histo_sort)
                fmo.seek(0)
                fmo.writelines(lines)

            print("check histos")
            exit(0)

            # run R script
            cmd = f"Rscript bin/SegStats.R AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerEqual{ak}.hist AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerOver{ak}.hist AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_AllMarkers.hist AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerSeg.png"
            subprocess.run(cmd, shell=True)

            print('check plots and histos')
            sys.exit(0)

            # perform analysis on progeny of G
            with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as ff1:
                for f1_prog in ff1:
                    f1_prog = f1_prog.strip().split()
                    if G in {f1_prog[3], f1_prog[4]}:
                        # check if call for F1 progeny exists
                        if not os.path.exists(f"AFLAP_tmp/04/Call/{f1_prog[0]}_{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.txt"):
                            raise FileNotFoundError(f"Error in 05_ObtainSegStats.py: Count for {f1_prog[0]} could not be found. Rerun 04_Genotyping.py.")
                            sys.exit(1)

                        # initialize variables
                        m_count = 0
                        cov = 0

                        # determine m_count based on F1 progeny's call file
                        with open(f"AFLAP_tmp/04/Call/{f1_prog[0]}_{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.txt", 'r') as fpcall:
                            for call in fpcall:
                                m_count += int(call.strip())

                        # determine cov based on F1 progeny's count file
                        with open(f"AFLAP_tmp/04/Count/{f1_prog[0]}_{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.txt", 'r') as fpcount:
                            c_dict = dict()
                            for count in fpcount:
                                if not len(count.strip()): continue     # skip empty lines
                                cval = int(count.strip().split()[1])

                                if cval not in c_dict: c_dict[cval] = 1
                                else: c_dict[cval] += 1

                            # find most frequent count
                            max = -np.inf
                            for c in c_dict:
                                if c_dict[c] > max:
                                    max = c_dict[c]
                                    cov = c

                            # confirm if cov's peak is 1
                            if cov == 1:
                                not1 = 0
                                is1 = 0
                                for c in c_dict:
                                    if c == 1: is1 = c_dict[c]
                                    if c > 5:  not1 += c_dict[c]

                                # remove data for coverage counts 0-5
                                if not1 > is1:
                                    c_dict.pop(1, 2, 3, 4, 5)

                                    # find most frequent count again
                                    max = -np.inf
                                    for c in c_dict:
                                        if c_dict[c] > max:
                                            max = c_dict[c]
                                            cov = c

                        print(f"{f1_prog[0]}\t{m_count}\t{cov}")




if __name__ == "__main__":
    main()