import argparse
import math
import os
import pandas as pd

from get_LA_info import get_LA_info
from seg_stats import get_seg_stats
from kmercov_x_markercount import plot_cov_and_mcount

#################################################
#       A shell script to obtain segregation statistics and exclude progeny which have low coverage.
#################################################

def histo_sort(line:str)->float:
    line_fields = line.strip().split()
    return float(line_fields[0])

def get_count_frequency(df:pd.DataFrame)->pd.DataFrame:
    return df.groupby("Frequency")["Frequency"].count().rename("Frequency Count").to_frame().reset_index(drop=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ObtainSegStats', description='A script to plot marker distributions in progeny.')
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--SDL', type=float, default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--SDU', type=float, default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    args = parser.parse_args()

    # initialize variables
    ak = 2 * int(args.kmer) - 1

    # analyze all parents whom we can create a genetic map for
    try:
        print("Performing segment statistics analysis...")
        list_of_Gs = get_LA_info()
        for G_info in list_of_Gs:
            G, LO, UP, P0 = G_info

            # check for num progs of G
            num_progs = 0
            with open("AFLAP_tmp/Crosses.txt", 'r') as fcrosses:
                for cross in fcrosses:
                    cross = cross.strip().split()

                    if G in (cross[2], cross[3]):
                        num_progs += int(cross[0])
            if not num_progs: raise ValueError("Invalid number of progeny.")

            print(f"\t\t{num_progs} Genotype calls for {G} detected. Summarizing...")

            # get marker stats
            tsv = pd.read_csv(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t')

            ## get frequencies (sum of calls / number of progeny)
            tsv["Frequency"] = tsv.iloc[:, 3:].sum(axis=1).div(num_progs)

            ## MarkerAll
            mal = get_count_frequency(tsv)
            ## MarkerEquals
            meq = tsv.loc[tsv["MarkerLength"].astype(int) == 61]
            meq = get_count_frequency(meq)
            ## MarkerOver
            mov = tsv.loc[tsv["MarkerLength"].astype(int) > 61]
            mov = get_count_frequency(mov)

            # get segment statistics
            get_seg_stats(mal, meq, mov, ak, f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}_MarkerSeg.png")

            # initialize dataframe of marker counts
            mc_df = pd.DataFrame(columns=["F1 Prog", "Marker Count", "K-mer Coverage"])

            # perform analysis on progeny of G
            with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as ff1:
                f1_prog_set = set()
                for f1_prog in ff1:
                    f1_prog = f1_prog.strip().split()
                    # skip over if f1_prog encountered already
                    if f1_prog[0] in f1_prog_set: continue
                    f1_prog_set.add(f1_prog[0])

                    if G in {f1_prog[3], f1_prog[4]}:
                        # check if call for F1 progeny exists
                        if not os.path.exists(f"AFLAP_tmp/04/Call/{f1_prog[0]}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt"):
                            raise FileNotFoundError(f"Count for {f1_prog[0]} could not be found. Rerun 04_Genotyping.py.")

                        # initialize variables
                        m_count = 0
                        cov = 0

                        # determine m_count based on F1 progeny's call file
                        with open(f"AFLAP_tmp/04/Call/{f1_prog[0]}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt", 'r') as fpcall:
                            for call in fpcall:
                                m_count += int(call.strip())

                        # determine cov based on F1 progeny's count file
                        with open(f"AFLAP_tmp/04/Count/{f1_prog[0]}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt", 'r') as fpcount:
                            c_dict = dict()
                            for count in fpcount:
                                if not len(count.strip()): continue     # skip empty lines
                                cval = int(count.strip().split()[1])

                                if cval not in c_dict: c_dict[cval] = 1
                                else: c_dict[cval] += 1

                            # delete 0 from c_dict
                            del c_dict[0]

                            # find most frequent count
                            max = -math.inf
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

                                # remove data for coverage counts 1 - 5
                                if not1 > is1:
                                    for i in range(1, 6):
                                        del c_dict[i]

                                    # find most frequent count again
                                    max = -math.inf
                                    for c in c_dict:
                                        if c_dict[c] > max:
                                            max = c_dict[c]
                                            cov = c

                        mc_df.loc[len(mc_df.index)] = [f1_prog[0], m_count, cov]

            # plot k-mer coverage and marker count
            plot_cov_and_mcount(mc_df, f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}_KmerCovXMarkerCount.png")

            if (not args.LOD == 2):
                print("\t\t\tAFLAP ran in low coverage mode. Coverage cut-off not run. Please manually remove any isolates you wish to exclude from $Ped and rerun AFLAP.\n" +
                    f"\t\t\tIt is possible that two peaks will be shown in AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}_MarkerSeg.png.\n" +
                    "\t\t\tIf that is the case please rerun AFLAP.sh providing -d and -D for lower and upper limits for marker filtering.")

            # filter out progeny with coverage < LOD
            low_cov = mc_df.loc[mc_df["K-mer Coverage"].astype(int) < int(args.LOD)]
            for i in low_cov.index:
                print(f"\t\t\t{low_cov['F1 Prog'][i]} appears to be low coverage. Will be excluded.")

            # create filtered tsv file
            tsv_filtered = tsv.loc[tsv["Frequency"].astype(float).between(args.SDL, args.SDU)]
            # remove Frequency column from tsv file
            tsv_filtered = tsv_filtered.iloc[:, :-1]
            tsv_filtered.to_csv(f"AFLAP_tmp/05/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv", sep='\t', index=False)

            print(f"\tFinished obtaining segment statistics for {G}.")
    
    except Exception as e:
        print(f"Error in 05_ObtainSegStats.py: {e}")
        exit(1)