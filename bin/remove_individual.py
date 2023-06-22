import pandas as pd

def remove_individual(rem_ind, pedigree)->None:
    # identify the individual to be removed
    ped_df = pd.read_csv(pedigree, sep='\t', header=None, usecols=[0,1], names=["Individual", "Generation"])
    if rem_ind not in ped_df["Individual"].unique():
        raise ValueError(f"Unable to remove {rem_ind} since it does not exist in the pedigree file.")

    ind_df = ped_df.loc[ped_df["Individual"] == rem_ind]
    # detect if more than one generation given to removed individual (FIXME: could have been identified in ped_analysis.py)
    if (ind_df["Generation"].unique() != 1):
        raise ValueError(f"{rem_ind} found to have been identified with more than one generation type.")
    rem_ind_gen = int(ind_df["Generation"].unique()[0])

    # remove individual from all analysis files
    ## figure out if we need to remove the individual's files from 01
    if rem_ind_gen == 0:
        # TODO: determine whether or not it is allowed to delete a parent
        # remove from 02
        # remove from 
        pass
    else:
        # TODO: remove from 01?
        pass
