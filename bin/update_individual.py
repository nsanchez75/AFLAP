import glob
import pandas as pd
import os

def remove_files(paths:list[str], parent_path:str)->None:
    tmpfiles = list()
    for path in paths:
        tmpfiles += glob.glob(path)
    if not tmpfiles:
        print(f"WARNING: No files have been removed from {parent_path}.")
    for tmpf in tmpfiles:
        os.remove(tmpf)

def update_individual(rem_ind, pedigree)->None:
    # identify the individual to be removed
    ped_df = pd.read_csv(pedigree, sep='\t', header=None, usecols=[0,1], names=["Individual", "Generation"])
    if rem_ind not in ped_df["Individual"].unique():
        exit(f"An error occurred: Unable to remove {rem_ind} since it does not exist in the pedigree file. Make sure you are using the same pedigree file from the previous run.")

    ind_df = ped_df.loc[ped_df["Individual"] == rem_ind]
    # detect if more than one generation given to removed individual (FIXME: could have been identified in ped_analysis.py)
    if (ind_df["Generation"].unique() != 1):
        raise ValueError(f"{rem_ind} found to have been identified with more than one generation type.")
    rem_ind_gen = int(ind_df["Generation"].unique()[0])

    # remove individual from all analysis files
    match rem_ind_gen:
        case 0:
            # remove jellyfish count
            remove_files([f"AFLAP_tmp/01/F0Count/{rem_ind}*"], "AFLAP_tmp/01")
            # remove histoplots
            remove_files([f"AFLAP_tmp/02/{rem_ind}*", f"AFLAP_tmp/02/F0Histo/{rem_ind}*"], "AFLAP_tmp/02")
            # remove markers
            remove_files([f"AFLAP_tmp/03/{rem_ind}*", f"AFLAP_tmp/03/F0Markers/{rem_ind}*"], "AFLAP_tmp/03")
            # remove tsv file
            remove_files([f"AFLAP_tmp/04/{rem_ind}*"], "AFLAP_tmp/04")
        case 1:
            # remove jellyfish count
            remove_files([f"AFLAP_tmp/01/F1Count/{rem_ind}*"], "AFLAP_tmp/01")
            # remove count and call
            remove_files([f"AFLAP_tmp/04/Call/{rem_ind}*", f"AFLAP_tmp/04/Count/{rem_ind}*"], "AFLAP_tmp/04")
        case 2:
            # remove jellyfish count
            remove_files([f"AFLAP_tmp/01/F2Count/{rem_ind}*"], "AFLAP_tmp/01")
            # remove count and call
            remove_files([f"AFLAP_tmp/04/Call/{rem_ind}*", f"AFLAP_tmp/04/Count/{rem_ind}*"], "AFLAP_tmp/04")