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
    print(ped_df)
    if not ped_df.loc[~ped_df["Generation"].astype(int).isin([0, 1, 2])].empty:
        exit("An error occurred: Individuals found with invalid generation type.")
    parents = ped_df.loc[ped_df["Generation"].astype(int) == 0]
    f1_progs = ped_df.loc[ped_df["Generation"].astype(int) == 1]
    f2_progs = ped_df.loc[ped_df["Generation"].astype(int) == 2]

    # order pedigrees by individual name
    parents = parents.sort_values(by="Individual")
    f1_progs = f1_progs.sort_values(by="Individual")
    f2_progs = f2_progs.sort_values(by="Individual")

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
    parents.to_csv("AFLAP_tmp/Pedigree_F0.txt", sep='\t', header=None, index=False)
    f1_progs.to_csv("AFLAP_tmp/Pedigree_F1.txt", sep='\t', header=None, index=False)
    f2_progs.to_csv("AFLAP_tmp/Pedigree_F2.txt", sep='\t', header=None, index=False)

    # perform LA analysis
    checked_parents = set()
    with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as fin, \
            open("AFLAP_tmp/LA.txt", 'w') as fla,       \
            open("AFLAP_tmp/noLA.txt", 'w') as fnola:

        for line in fin:
            line = line.strip().split()
            if line[0] not in checked_parents:
                if "NA" in (line[3], line[4]):
                    fnola.write(f"{line[0]}\n")
                    checked_parents.add(line[0])
                elif [isinstance(x, int) for x in [line[3], line[4]]]:
                    if int(line[3]) > int(line[4]):
                        print(f"{line[3]} | {line[4]}")
                        exit("An error occurred: Cannot have a lower bound higher than an upper bound.")

                    # write {parent} {lower bound} {upper bound} to LA.txt
                    fla.write(f"{line[0]}\t{line[3]}\t{line[4]}\n")
                    checked_parents.add(line[0])
                else:
                    exit("An error occurred: Invalid bound entry in Pedigree_F0.txt.")

    print("Finished pedigree file analysis.")
