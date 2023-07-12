import pandas as pd
import re
import shutil

def write_cross(prog_info:pd.DataFrame, ftype:int, parents:list)->None:
    # check if parents are valid
    if (prog_info["MP"].astype(str) == "NA").any() | (prog_info["FP"].astype(str) == "NA").any():
        exit(f"An error occurred: There is an F{ftype} progeny that has 'NA' parent(s).")
    if (~prog_info["MP"].astype(str).isin(parents)).any() or (~prog_info["FP"].astype(str).isin(parents)).any():
        exit(f"An error occurred: An F{ftype} progeny comes from parent not found in the pedigree file.")
    if (prog_info["MP"].astype(str) == prog_info["FP"].astype(str)).any():
        exit(f"An error occurred: Identical crossed parents for an F{ftype} progeny identified.")

    # return cross info
    count_crosses = prog_info.groupby(["MP", "FP"]).size().reset_index()
    crosses = list(tuple(row) for row in count_crosses.itertuples(index=False))
    with open("AFLAP_tmp/Crosses.txt", 'r+') as fc:
         for cross in crosses:
              MP, FP, COUNT = cross
              fc.write(f"{COUNT}\t{ftype}\t{MP}\t{FP}\n")

def pedigree_analysis(pedigree: str)->None:
    # copy pedigree file into AFLAP_Results
    shutil.copy2(pedigree, "AFLAP_Results/Pedigree.txt")

    # remove comments in pedigree files (start with '#)
    with open(pedigree, 'r+') as f:
        lines = f.readlines()
        lines = [line for line in lines if not re.match(r'^\s*#', line)]
        f.seek(0)
        f.writelines(lines)

    # store pedigree filename to pedigree info file
    with open("AFLAP_tmp/PedigreeInfo.txt", 'w') as fpinfo:
        fpinfo.write(f"Source: {pedigree}")

    # create pedigree dataframe and categorize by generations
    ped_df = pd.read_csv(pedigree, sep='\t', header=None, names=["Individual", "Generation", "Path", "MP", "FP"])
    if not ped_df.loc[~ped_df["Generation"].astype(int).isin([0, 1, 2])].empty:
        exit("An error occurred: Individuals found with invalid generation type.")
    parents = ped_df.loc[ped_df["Generation"].astype(int) == 0].rename(columns={"MP": "LB", "FP": "UB"})
    f1_progs = ped_df.loc[ped_df["Generation"].astype(int) == 1]
    f2_progs = ped_df.loc[ped_df["Generation"].astype(int) == 2]

    # order pedigrees by individual name
    parents = parents.sort_values(by="Individual")
    f1_progs = f1_progs.sort_values(by="Individual")
    f2_progs = f2_progs.sort_values(by="Individual")

    # check if same individual does not have different parents
    for df in [f1_progs, f2_progs]:
        for ind in df["Individual"].unique():
            if df[df["Individual"] == ind]["MP"].unique().size != 1 or \
               df[df["Individual"] == ind]["FP"].unique().size != 1:
                exit("An error occurred: More than one parent detected for an individual.")

    # initialize variables
    p_list = parents["Individual"].unique()

    # count number of parents
    p_len = len(p_list)
    if p_len < 2:
        exit("An error occurred: Pedigree file does not have at least 2 parents.")
    elif p_len == 2: print("2 parents detected. This will be easy! Identifying cross(es)...")
    else: print(f"{p_len} parents detected. This will not be so easy! Identifying crosses...")

    # check prog pedigrees and write to cross file
    with open("AFLAP_tmp/Crosses.txt", 'w') as fc:
        pass
    write_cross(f1_progs, 1, p_list)
    write_cross(f2_progs, 2, p_list)

    # create categorized pedigree files
    parents.to_csv("AFLAP_tmp/Pedigree_F0.txt", sep='\t', header=True, index=False)
    f1_progs.to_csv("AFLAP_tmp/Pedigree_F1.txt", sep='\t', header=True, index=False)
    f2_progs.to_csv("AFLAP_tmp/Pedigree_F2.txt", sep='\t', header=True, index=False)

    # perform LA analysis
    with open("AFLAP_tmp/LA.txt", 'w') as fla, open("AFLAP_tmp/noLA.txt", 'w') as fnola:
        for parent in parents["Individual"].unique():
            parent_df = parents[parents["Individual"] == parent]
            if parent_df["LB"].unique().size != 1:
                exit(f"An error occurred: There must be only one lower bound for {parent}.")
            if parent_df["UB"].unique().size != 1:
                exit(f"An error occurred: There must be only one upper bound for {parent}.")

            lower_bound = parent_df["LB"].unique()[0]
            upper_bound = parent_df["UB"].unique()[0]
            if "NA" in (lower_bound, upper_bound):
                fnola.write(f"{parent}\n")
            elif [isinstance(x, int) for x in [lower_bound, upper_bound]]:
                if int(lower_bound) > int(upper_bound):
                    exit("An error occurred: Cannot have a lower bound higher than an upper bound.")
                fla.write(f"{parent}\t{lower_bound}\t{upper_bound}\n")
            else:
                exit("An error occurred: Invalid bound entry in Pedigree_F0.txt.")

    print("Finished pedigree file analysis.")
