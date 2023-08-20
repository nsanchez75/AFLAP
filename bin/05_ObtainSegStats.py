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
    return df.groupby("Frequency")["Frequency"].count() \
             .rename("Frequency Count").to_frame().reset_index(drop=False)

def progeny_analysis(progs_df:pd.DataFrame, f_type:str, G_info:tuple, mc_df:pd.DataFrame, kmer:int)->pd.DataFrame:
    G, LO, UP, P0, SEX = G_info

    for prog in progs_df["Individual"].unique().tolist():
        prog_df = progs_df[progs_df["Individual"].astype(str) == prog]
        if prog_df[(prog_df["MP"].astype(str) == G) | (prog_df["MP"].astype(str) == G)].empty: continue
        prog_call_file = f"AFLAP_tmp/04/{f_type}/Call/{prog}_{G}_m{kmer}_L{LO}_U{UP}_{P0}.txt"
        if not os.path.exists(prog_call_file):
            exit(f"An error occurred: Count for {prog} could not be found. Rerun 04_Genotyping.py.")

        # find progeny's marker count
        with open(prog_call_file, 'r') as fcall:
            call_vals = [int(line.strip()) for line in fcall]
            marker_count = sum(call_vals)

        # find progeny's coverage value
        coverage = 0
        count_df = pd.read_csv(f"AFLAP_tmp/04/{f_type}/Count/{prog}_{G}_m{kmer}_L{LO}_U{UP}_{P0}.txt", sep=' ', header=None, names=["Sequence", "Count"])
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

def filter_f1(kmer:int)->None:
    for G_info in get_LA_info():
        G, LO, UP, P0, SEX = G_info

        # check number of calls for G
        call_files = glob.glob("AFLAP_tmp/04/F1/Call/*.txt")
        call_files = list(filter(lambda x: True if (x[21:].split('_')[1] == G) else False, call_files))
        num_progs = len(call_files)
        if not num_progs: exit("An error occurred: Invalid number of progeny.")
        print(f"{num_progs} Genotype calls for {G} detected. Summarizing...")

        # perform analysis on F1 progeny of G
        mc_df = pd.DataFrame(columns=["Prog", "Marker Count", "K-mer Coverage"])
        f1_progs_df = pd.read_csv("AFLAP_tmp/Pedigree_F1.txt", sep='\t')
        mc_df = progeny_analysis(f1_progs_df, "F1", G_info, mc_df, kmer)

        # plot k-mer coverage and marker count
        plot_cov_and_mcount(mc_df, f"AFLAP_Results/{G}_m{kmer}_L{LO}_U{UP}_{P0}_KmerCovXMarkerCount.png")

        # filter genotype tables
        genotype_file = f"AFLAP_tmp/04/{G}_F1_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv"
        if not os.path.exists(genotype_file):
            print(f"Genotype for F1 progeny of {G} not detected. Skipping.")
            continue

        tsv_df = pd.read_csv(f"AFLAP_tmp/04/{G}_F1_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t')
        tsv_df["Frequency"] = tsv_df.iloc[:, 3:].sum(axis=1).div(num_progs)

        marker_all = get_count_frequency(tsv_df)
        marker_equals = get_count_frequency(tsv_df[tsv_df["MarkerLength"].astype(int) == 61])
        marker_over = get_count_frequency(tsv_df[tsv_df["MarkerLength"].astype(int) > 61])

        # get marker statistics
        seg_png = f"AFLAP_Results/{G}_m{kmer}_L{LO}_U{UP}_{P0}_MarkerSeg.png"
        ak = 2 * kmer - 1
        get_seg_stats(marker_all, marker_equals, marker_over, ak, seg_png)

        if not args.LOD == 2:
            print("\tAFLAP ran in low coverage mode. Coverage cut-off not run. Please manually remove any isolates you wish to exclude from the pedigree file and rerun AFLAP.\n" +
                 f"\tIt is possible that two peaks will be shown in AFLAP_Results/{G}_m{kmer}_L{LO}_U{UP}_{P0}_MarkerSeg.png.\n" +
                  "\tIf that is the case please rerun AFLAP.sh providing -d and -D for lower and upper limits for marker filtering.")

        # filter out progeny with coverage < LOD
        low_cov = mc_df[mc_df["K-mer Coverage"].astype(int) < int(args.LOD)]
        for i in low_cov.index:
            # TODO: exclude prog from tsv
            print(f"\t{low_cov['Prog'][i]} appears to be low coverage. Will be excluded.")

        # create filtered tsv file
        tsv_df_filtered = tsv_df[tsv_df["Frequency"].astype(float).between(args.SDL, args.SDU)]
        ## remove Frequency column from tsv file
        tsv_df_filtered = tsv_df_filtered.iloc[:, :-1]
        tsv_df_filtered.to_csv(f"AFLAP_tmp/05/{G}_F1_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv", sep='\t', index=False)

        print(f"Finished obtaining F1 segment statistics for {G}.")

def agg_function(series:pd.Series):
    unique_values = series.astype(str).unique()
    if len(unique_values) == 1:
        return unique_values[0]
    else:
        return ''.join(unique_values)

def filter_f2(kmer:int, xx_filter:float)->None:
    # get male and female marker dataframes
    parents = ['', '']
    for G_info in get_LA_info():
        G, LO, UP, P0, SEX = G_info
        marker_df = pd.read_csv(f"AFLAP_tmp/04/{G}_F2_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t').drop(columns=["MarkerLength"])

        if SEX == "male":
            male_marker_df = marker_df
            parents[0] = G
        else:
            female_marker_df = marker_df
            parents[1] = G

    if male_marker_df.empty or female_marker_df.empty:
        exit("An error occurred: There should be two marker tables for the female and male parent. To fix this, rerun 04_Genotyping.py.")

    # combine dataframes to make filtered F2 dataframe
    comb_marker_df = pd.concat([male_marker_df, female_marker_df], ignore_index=True)
    comb_marker_df = comb_marker_df.groupby(["MarkerSequence", "MarkerID"], as_index=False).agg({col: agg_function for col in comb_marker_df.columns[2:]})
    replace_dict = {'A': 'AA', 'B': 'BB', 'X': 'XX', 'AX': 'AA', 'XB': 'BB'}
    comb_marker_df[comb_marker_df.columns[2:].tolist()] = comb_marker_df[comb_marker_df.columns[2:].tolist()].replace(replace_dict)

    # get sequence statistics
    freq_cols = list()
    freq_df = pd.DataFrame()
    for val in ['XX', 'AA', 'BB', 'AB']:
        freq_cols.append(pd.Series(comb_marker_df.iloc[:, 2:].eq(val).sum(axis=1) / (comb_marker_df.shape[1] - 2), name=f"{val} Frequency"))
    freq_df = pd.concat(freq_cols, axis=1)
    comb_marker_df = pd.concat([comb_marker_df, freq_df], axis=1)

    # create filtered table
    comb_marker_df.to_csv(f"AFLAP_tmp/05/{parents[0]}x{parents[1]}_F2_m{kmer}.Genotypes.MarkerID.Filtered.tsv", sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ObtainSegStats', description='A script to plot marker distributions in progeny.')
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--SDL', type=float, default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--SDU', type=float, default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    parser.add_argument('-f', '--fXX', type=float, default=0.1, help='Limit for how many XX can exist in a row. If surpassed then sequence is not considered for analysis. Default [0.1].')
    args = parser.parse_args()

    # create directory
    os.makedirs("AFLAP_tmp/05", exist_ok=True)

    # analyze all parents whom we can create a genetic map for
    print("Performing segment statistics analysis...")
    list_of_Gs = get_LA_info()
    for f_type in ["F1", "F2"]:
        if not len(os.listdir(f"AFLAP_tmp/01/{f_type}Count")):
            print(f"No {f_type} progeny found. Skipping.")
            continue

        if f_type == "F1": filter_f1(args.kmer)
        else:              filter_f2(args.kmer, args.fXX)
