import argparse
import glob
import math
import os
import pandas as pd

from get_LA_info import get_LA_info
from seg_stats import get_seg_stats
from kmercov_x_markercount import plot_cov_and_mcount

#################################################
#       A shell script to obtain segregation statistics and exclude progeny which have low coverage.
#################################################

def get_count_frequency(df:pd.DataFrame)->pd.DataFrame:
    return df.groupby("Frequency")["Frequency"].count().rename("Frequency Count").to_frame().reset_index(drop=False)

def progeny_analysis(progs_df:pd.DataFrame, f_type:str, G_info:tuple, mc_df:pd.DataFrame)->pd.DataFrame:
    G, LO, UP, P0, SEX = G_info

    for prog in progs_df["Individual"].unique().tolist():
        prog_df = progs_df[progs_df["Individual"].astype(str) == prog]
        if prog_df[(prog_df["MP"].astype(str) == G) | (prog_df["MP"].astype(str) == G)].empty: continue
        prog_call_file = f"AFLAP_tmp/04/{f_type}/Call/{prog}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt"
        if not os.path.exists(prog_call_file):
            exit(f"An error occurred: Count for {prog} could not be found. Rerun 04_Genotyping.py.")

        # find progeny's marker count
        with open(prog_call_file, 'r') as fcall:
            call_vals = [int(line.strip()) for line in fcall]
            marker_count = sum(call_vals)

        # find progeny's coverage value
        coverage = 0
        count_df = pd.read_csv(f"AFLAP_tmp/04/{f_type}/Count/{prog}_{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.txt", sep=' ', header=None, names=["Sequence", "Count"])
        c_dict = count_df.groupby("Count").count().to_dict()["Sequence"]
        del c_dict[0]   # remove 0 from dictionary

        ## find most frequent count
        max = -math.inf
        for c in c_dict:
            if c_dict[c] > max:
                max = c_dict[c]
                coverage = c

        ## confirm if coverage peak is 1
        if coverage == 1:
            not1 = is1 = 0
            for i in range(2, 6): del c_dict[i]
            for c in c_dict:
                if c == 1: is1 = c_dict[c]
                not1 += c_dict[c]

            if not1 > is1:
                del c_dict[1]
                # find most frequent count again
                max = -math.inf
                for c in c_dict:
                    if c_dict[c] > max:
                        max = c_dict[c]
                        coverage = c

        mc_df.loc[len(mc_df.index)] = [prog, marker_count, coverage]

    return mc_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ObtainSegStats', description='A script to plot marker distributions in progeny.')
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--SDL', type=float, default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--SDU', type=float, default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    args = parser.parse_args()

    # create directory
    os.makedirs("AFLAP_tmp/05", exist_ok=True)

    # analyze all parents whom we can create a genetic map for
    print("Performing segment statistics analysis...")
    list_of_Gs = get_LA_info()
    for G_info in list_of_Gs:
        G, LO, UP, P0, SEX = G_info

        # check number of calls for G
        call_files = glob.glob("AFLAP_tmp/04/F*/Call/*.txt")
        call_files = [file[21:] for file in call_files]
        call_files = list(filter(lambda x: True if (x.split('_')[1] == G) else False, call_files))
        num_progs = len(call_files)
        if not num_progs: exit("An error occurred: Invalid number of progeny.")
        print(f"{num_progs} Genotype calls for {G} detected. Summarizing...")

        # perform analysis on progeny of G
        mc_df = pd.DataFrame(columns=["Prog", "Marker Count", "K-mer Coverage"])
        f1_progs_df = pd.read_csv("AFLAP_tmp/Pedigree_F1.txt", sep='\t')
        mc_df = progeny_analysis(f1_progs_df, "F1", G_info, mc_df)
        f2_progs_df = pd.read_csv("AFLAP_tmp/Pedigree_F2.txt", sep='\t')
        mc_df = progeny_analysis(f2_progs_df, "F2", G_info, mc_df)

        # plot k-mer coverage and marker count
        plot_cov_and_mcount(mc_df, f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}_KmerCovXMarkerCount.png")

        # filter genotype tables
        for f_type in ["F1", "F2"]:
            genotype_file = f"AFLAP_tmp/04/{G}_{f_type}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv"
            if not os.path.exists(genotype_file):
                print(f"Genotype for {f_type} progeny of {G} not detected. Skipping.")
                continue

            tsv_df = pd.read_csv(f"AFLAP_tmp/04/{G}_{f_type}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t')
            tsv_df["Frequency"] = tsv_df.iloc[:, 3:].sum(axis=1).div(num_progs)

            marker_all = get_count_frequency(tsv_df)
            marker_equals = get_count_frequency(tsv_df[tsv_df["MarkerLength"].astype(int) == 61])
            marker_over = get_count_frequency(tsv_df[tsv_df["MarkerLength"].astype(int) > 61])

            # get marker statistics
            seg_png = f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}_MarkerSeg.png"
            ak = 2 * int(args.kmer) - 1
            get_seg_stats(marker_all, marker_equals, marker_over, ak, seg_png)


            if not args.LOD == 2:
                print("\tAFLAP ran in low coverage mode. Coverage cut-off not run. Please manually remove any isolates you wish to exclude from $Ped and rerun AFLAP.\n" +
                     f"\tIt is possible that two peaks will be shown in AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}_MarkerSeg.png.\n" +
                      "\tIf that is the case please rerun AFLAP.sh providing -d and -D for lower and upper limits for marker filtering.")

            # filter out progeny with coverage < LOD
            low_cov = mc_df[mc_df["K-mer Coverage"].astype(int) < int(args.LOD)]
            for i in low_cov.index:
                print(f"\t{low_cov['Prog'][i]} appears to be low coverage. Will be excluded.")

            # create filtered tsv file
            tsv_df_filtered = tsv_df[tsv_df["Frequency"].astype(float).between(args.SDL, args.SDU)]
            # remove Frequency column from tsv file
            tsv_df_filtered = tsv_df_filtered.iloc[:, :-1]

            # identify identical loci sequences if working on F2
            if f_type == "F2":
                set_of_loci_seqs = set()
                with open("AFLAP_tmp/03/SimGroups/identical_loci.txt", 'r') as floci:
                    [set_of_loci_seqs.add(line.strip()) for line in floci]
                condition = tsv_df_filtered["MarkerSequence"].astype(str).isin(set_of_loci_seqs)
                replaced_columns = tsv_df_filtered.columns.tolist()[3:]
                tsv_df_filtered.loc[condition, replaced_columns] = 2

            tsv_df_filtered.to_csv(f"AFLAP_tmp/05/{G}_{f_type}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv", sep='\t', index=False)

            print(f"Finished obtaining {f_type} segment statistics for {G}.")
