import argparse
import glob
import math
import os
import pandas as pd

from get_LA_info import get_LA_info
from seg_stats import get_seg_stats
from kmercov_x_markercount import plot_cov_and_mcount

#################################################
#       A Python script to obtain segregation statistics and exclude progeny which have low coverage.
#################################################

def get_count_frequency(df:pd.DataFrame)->pd.DataFrame:
    return df.groupby("Frequency")["Frequency"].count().rename("Frequency Count").to_frame().reset_index(drop=False)

def progeny_analysis(progs_df:pd.DataFrame, f_type:str, G_info:tuple, mc_df:pd.DataFrame, kmer:int)->pd.DataFrame:
    G, LO, UP, P0, _ = G_info

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
        count_df = pd.read_csv(f"AFLAP_tmp/04/{f_type}/Count/{prog}_{G}_m{kmer}_L{LO}_U{UP}_{P0}.txt", sep=' ',
                               header=None, names=["Sequence", "Count"])
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

def genotype_table_stats(marker_df:pd.DataFrame, G_info:tuple, f_type:str, kmer:int, LOD:int, SDL:float, SDU:float)->tuple[pd.DataFrame, set]:
    G, LO, UP, P0, _ = G_info
    filtered_progs = set()

    # check number of calls for parent
    call_files = glob.glob(f"AFLAP_tmp/04/{f_type}/Call/*.txt")
    call_files = list(filter(lambda x: True if (x[21:].split('_')[1] == G) else False, call_files))
    num_progs = len(call_files)
    if not num_progs: exit("An error occurred: Invalid number of progeny.")
    print(f"{num_progs} Genotype calls for {G} detected. Summarizing...")

    # perform analysis on F2 progeny of G
    mc_df = pd.DataFrame(columns=["Prog", "Marker Count", "K-mer Coverage"])
    progs_df = pd.read_csv(f"AFLAP_tmp/Pedigree_{f_type}.txt", sep='\t')
    mc_df = progeny_analysis(progs_df, f_type, G_info, mc_df, kmer)

    # plot k-mer coverage and marker count
    plot_cov_and_mcount(mc_df, f"AFLAP_Results/{G}_{f_type}_m{kmer}_L{LO}_U{UP}_{P0}_KmerCovXMarkerCount.png")

    # get marker statistics
    if f_type == "F1": marker_df["Frequency"] = marker_df.iloc[:, 3:].astype(int).sum(axis=1).div(num_progs)
    else:              marker_df["Frequency"] = marker_df.iloc[:, 3:].astype(str).ne('X').sum(axis=1).div(num_progs)
    marker_all = get_count_frequency(marker_df)
    marker_equals = get_count_frequency(marker_df[marker_df["MarkerLength"].astype(int) == 61])
    marker_over = get_count_frequency(marker_df[marker_df["MarkerLength"].astype(int) > 61])
    seg_png = f"AFLAP_Results/{G}_{f_type}_m{kmer}_L{LO}_U{UP}_{P0}_MarkerSeg.png"
    ak = 2 * kmer - 1
    get_seg_stats(marker_all, marker_equals, marker_over, ak, seg_png)

    # filter out progeny with coverage < LOD
    if not LOD == 2:
        print("\tAFLAP ran in low coverage mode. Coverage cut-off not run. Please manually remove any isolates you wish to exclude from the pedigree file and rerun AFLAP.\n" +
            f"\tIt is possible that two peaks will be shown in {seg_png}.\n" +
            "\tIf that is the case please rerun AFLAP.py providing -d and -D for lower and upper limits for marker filtering.")
    low_cov = mc_df[mc_df["K-mer Coverage"].astype(int) < LOD]
    for i in low_cov.index:
        filtered_progs.add(low_cov['Prog'][i])
        print(f"\t{low_cov['Prog'][i]} appears to be low coverage. Will be excluded.")

    # drop sequences outside of frequency bounds
    marker_df = marker_df[marker_df["Frequency"].astype(float).between(SDL, SDU)]
    marker_df = marker_df.drop(columns=["Frequency"])

    return (marker_df, filtered_progs)

def filter_f1(kmer:int, LOD:int, SDL:float, SDU:float)->None:
    for G_info in get_LA_info():
        G, LO, UP, P0, _ = G_info
        marker_df = pd.read_csv(f"AFLAP_tmp/04/{G}_F1_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t')

        # get genotype table statistics
        marker_df, filtered_progs = genotype_table_stats(marker_df, G_info, "F1", kmer, LOD, SDL, SDU)

        # remove LOD-filtered progeny from male and female dataframes
        marker_df = marker_df.drop(columns=list(filtered_progs))

        # combine marker ID and marker length
        marker_df["MarkerID"] = marker_df["MarkerID"].astype(str) + '_' + marker_df["MarkerLength"].astype(str)
        marker_df = marker_df.drop(columns=["MarkerLength"])

        marker_df.to_csv(f"AFLAP_tmp/05/{G}_F1_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv", sep='\t', index=False)
        print(f"Finished creating F1 genotype table for {G}.")

def agg_function(series:pd.Series)->str:
    unique_values = series.astype(str).unique()
    return unique_values[0] if (len(unique_values)) else ''.join(unique_values)

def filter_f2(kmer:int, LOD:int, SDL:float, SDU:float, xx_filter:float=None)->None:
    # get male and female marker dataframes
    parents = [[], []]
    filtered_progs = set()
    for G_info in get_LA_info():
        G, LO, UP, P0, SEX = G_info
        marker_df = pd.read_csv(f"AFLAP_tmp/04/{G}_F2_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.tsv", sep='\t')

        # get genotype table statistics
        marker_df, new_filtered_progs = genotype_table_stats(marker_df, G_info, "F2", kmer, LOD, SDL, SDU)
        filtered_progs.union(new_filtered_progs)

        # drop marker length (not used in F2's filtered table)
        marker_df = marker_df.drop(columns=["MarkerLength"])

        # extract necessary info for filtered table
        if SEX == "male":
            male_marker_df = marker_df
            parents[0].append(G)
        else:
            female_marker_df = marker_df
            parents[1].append(G)

    # check if dataframes contain any sequences
    if male_marker_df.empty or female_marker_df.empty:
        exit("An error occurred: There should be two marker tables for the female and male parent. To fix this, rerun 04_Genotyping.py.")

    # concatenate all parents of same sex together
    parents = ['_'.join(p) for p in parents]

    # remove LOD-filtered progeny from male and female dataframes
    for df in [male_marker_df, female_marker_df]: df = df.drop(columns=list(filtered_progs))

    # combine dataframes to make filtered F2 dataframe
    comb_marker_df = pd.concat([male_marker_df, female_marker_df], ignore_index=True)
    comb_marker_df = comb_marker_df.groupby(["MarkerSequence", "MarkerID"], as_index=False).agg({col: agg_function for col in comb_marker_df.columns[2:]})
    replace_dict = {'A': 'AA', 'B': 'BB', 'X': 'XX', 'AX': 'AA', 'XB': 'BB'}
    comb_marker_df[comb_marker_df.columns[2:].tolist()] = comb_marker_df[comb_marker_df.columns[2:].tolist()].replace(replace_dict)

    # get sequence statistics
    freq_df = pd.DataFrame()
    freq_cols = list()
    for val in ['XX', 'AA', 'BB', 'AB']:
        freq_cols.append(pd.Series(comb_marker_df.iloc[:, 2:].astype(str).eq(val).sum(axis=1).div(comb_marker_df.shape[1] - 2), name=f"{val} Frequency"))
    freq_df = pd.concat(freq_cols, axis=1)
    comb_marker_df = pd.concat([comb_marker_df, freq_df], axis=1)

    # if necessary, filter out sequences by XX frequency
    if xx_filter is not None: comb_marker_df = comb_marker_df[comb_marker_df["XX Frequency"].astype(float) <= xx_filter]
    else:
        print("XX Filter not specified. Using median of all XX Frequency values as a filter...")
        comb_marker_df = comb_marker_df[comb_marker_df["XX Frequency"].astype(float) <= comb_marker_df["XX Frequency"].astype(float).quantile(0.75)]

    # create frequency stats table
    freq_stats_df = comb_marker_df[["MarkerSequence", "MarkerID", "XX Frequency", "AA Frequency", "BB Frequency", "AB Frequency"]]
    freq_stats_df.to_csv(f"AFLAP_tmp/05/{parents[0]}x{parents[1]}_F2_m{kmer}.FilteredFrequencyStats.tsv", sep='\t', index=False)
    comb_marker_df = comb_marker_df.drop(columns=["XX Frequency", "AA Frequency", "BB Frequency", "AB Frequency"])

    # create filtered genotype table
    comb_marker_df.to_csv(f"AFLAP_tmp/05/{parents[0]}x{parents[1]}_F2_m{kmer}.Genotypes.MarkerID.Filtered.tsv", sep='\t', index=False)
    print("Finished creating F2 filtered genotype table.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ObtainSegStats', description='A script to plot marker distributions in progeny.')
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--SDL', type=float, default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--SDU', type=float, default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    parser.add_argument('-f', '--fXX', type=float, default=None, help='Limit for how many XX can exist in a row. If surpassed then sequence is not considered for analysis. Default [None].')
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

        if f_type == "F1": filter_f1(args.kmer, args.LOD, args.SDL, args.SDU)
        else:              filter_f2(args.kmer, args.LOD, args.SDL, args.SDU, args.fXX)
