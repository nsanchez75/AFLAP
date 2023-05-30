import pandas as pd
import shutil

def check_prog(prog_info:list, parents:set, progs:set, cross_dict:dict)->None:
    # parents unidentified
    if prog_info[3] not in parents or prog_info[4] not in parents:
        raise ValueError(f"Progeny {prog_info[0]} comes from parent not found in the pedigree file.")
    # unknown parent
    if "NA" in (prog_info[3], prog_info[4]):
        raise ValueError(f"Progeny {prog_info[0]} has 'NA' parent(s).")
    # crossing identical parent
    if prog_info[3] == prog_info[4]:
        raise ValueError(f"Identical crossed parents identified for {prog_info[0]}.")

    if prog_info[0] not in progs:
        progs.add(prog_info[0])
        # create new cross in crosses dictionary
        cross = f"{prog_info[3]} {prog_info[4]}"
        if cross not in cross_dict:
            if f"{prog_info[4]} {prog_info[3]}" in cross_dict:
                raise ValueError("A parent in the pedigree is being treated as both male and female.")
            cross_dict[cross] = 1
        else: cross_dict[cross] += 1

def sort_ped(f_type:str)->None:
    with open(f"AFLAP_tmp/Pedigree_{f_type}.txt", 'r+') as f:
        lines = sorted(f.readlines())
        f.seek(0)
        f.writelines(lines)

def pedigree_analysis(pedigree: str)->None:
    # try:
        # copy pedigree file into AFLAP_Results
        shutil.copy2(pedigree, "AFLAP_Results/Pedigree.txt")
        
        # create pedigree dataframe
        ped_df = pd.read_csv(pedigree, sep='\t', header=None, names=["Individual", "Generation", "Path", "MP", "FP"])

        # separate pedigree dataframe by generation
        if not ped_df.loc[~ped_df["Generation"].astype(int).isin([0, 1, 2])].empty:
            raise ValueError("Individuals found with invalid generation type.")
        parents = ped_df.loc[ped_df["Generation"].astype(int) == 0]
        f1_progs = ped_df.loc[ped_df["Generation"].astype(int) == 1]
        f2_progs = ped_df.loc[ped_df["Generation"].astype(int) == 2]

        # create LA.txt from parents
        la_df = parents.loc[(parents["MP"] != "NA") & (parents["FP"] != "NA")]
        nola_df = parents.loc[~((parents["MP"] != "NA") & (parents["FP"] != "NA"))]
        print("MP:")
        print(la_df.groupby("Individual")["MP"].nunique())
        print("FP:")
        print(la_df.groupby("Individual")["FP"].nunique())
        result = la_df.groupby("Individual").apply(lambda x : (x["MP"].nunique == 1) and (x["FP"].nunique == 1))
        print(result)
        # print(la_df)


        # create categorized pedigree files
        # parents.to_csv("AFLAP_tmp/Pedigree_F0.txt", sep='\t', header=None, index=False)
        # f1_progs.to_csv("AFLAP_tmp/Pedigree_F1.txt", sep='\t', header=None, index=False)
        # f2_progs.to_csv("AFLAP_tmp/Pedigree_F2.txt", sep='\t', header=None, index=False)

        # initialize cross dictionaries
        f1_crosses = dict()
        f2_crosses = dict()



        # try:
        #     # initialize variables
        #     parents    = set()
        #     f1_progs   = set()
        #     f2_progs   = set()
        #     f1_crosses = dict()
        #     f2_crosses = dict()

        #     # copy pedigree file into AFLAP_Results
        #     shutil.copy2(pedigree, "AFLAP_Results/Pedigree.txt")

        #     # FIXME: convert pedigree into pandas df to grab stuff faster
        #     # categorize pedigree into F0, F1, F2
        #     with open(pedigree, 'r') as fin,                        \
        #          open("AFLAP_tmp/PedigreeInfo.txt", 'w') as fped,   \
        #          open("AFLAP_tmp/Pedigree_F0.txt", 'w') as f0,      \
        #          open("AFLAP_tmp/Pedigree_F1.txt", 'w') as f1,      \
        #          open("AFLAP_tmp/Pedigree_F2.txt", 'w') as f2:

        #         # write pedigree file into Pedigree.txt
        #         fped.write(f"Source: {pedigree}")

        #         # categorize individuals by generation type
        #         for line in fin:
        #             cols = line.strip().split()
        #             # check line size
        #             if len(cols) != 5:
        #                 raise ValueError("Line in pedigree file over 5 columns.")

        #             if cols[1] == '0':
        #                 if cols[0] not in parents:
        #                     parents.add(cols[0])

        #                 f0.write(line)
        #             elif cols[1] == '1':
        #                 print(cols)
        #                 print(parents)
        #                 check_prog(cols, parents, f1_progs, f1_crosses)
        #                 f1.write(line)
        #             elif cols[1] == '2':
        #                 check_prog(cols, parents, f2_progs, f2_crosses)
        #                 f2.write(line)
        #             else:
        #                 raise ValueError("Pedigree file contains individual that is not F0, F1, or F2.")

        #     # sort pedigree files
        #     sort_ped("F0")
        #     sort_ped("F1")
        #     sort_ped("F2")

        #     # determine parent count
        #     if len(parents) < 2:
        #         raise ValueError("Pedigree file does not have at least 2 parents.")
        #     elif len(parents) == 2: print("2 parents detected. This will be easy! Identifying cross(es)...")
        #     else: print(f"{len(parents)} parents detected. This will not be so easy! Identifying crosses...")

        #     # identify and print cross
        #     with open("AFLAP_tmp/Crosses.txt", 'w') as f:
        #         ## F1 progeny
        #         print("F1 crosses that have been identified:")
        #         for cross in f1_crosses:
        #             c_vals = cross.strip().split()
        #             print(f"\t{c_vals[0]}x{c_vals[1]}")
        #             f.write(f"{f1_crosses[cross]} 1 {c_vals[0]} {c_vals[1]}")
        #         print()
        #         ## F2 progeny
        #         # TODO: implement stuff below when working with F2
        #         print("F2 crosses that have been identified:")
        #         for crosses in f2_crosses:
        #             c_vals = cross.strip().split()
        #             print(f"\t{c_vals[0]}x{c_vals[1]}")
        #             f.write(f"{f2_crosses[crosses]} 2 {c_vals[0]} {c_vals[1]}")

        #     # perform LA analysis
        #     checked_parents = set()
        #     with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as fin, \
        #          open("AFLAP_tmp/LA.txt", 'w') as fla,       \
        #          open("AFLAP_tmp/noLA.txt", 'w') as fnola:

        #         for line in fin:
        #             line = line.strip().split()
        #             if line[0] not in checked_parents:
        #                 if "NA" in (line[3], line[4]):
        #                     fnola.write(f"{line[0]}\n")
        #                     checked_parents.add(line[0])
        #                 elif [isinstance(x, int) for x in [line[3], line[4]]]:
        #                     if line[3] > line[4]:
        #                         raise("Cannot have a lower bound higher than an upper bound.")

        #                     # write {parent} {lower bound} {upper bound} to LA.txt
        #                     fla.write(f"{line[0]} {line[3]} {line[4]}\n")
        #                     checked_parents.add(line[0])
        #                 else:
        #                     raise ValueError("Invalid bound entry in Pedigree_F0.txt.")

        #     print("Finished pedigree file analysis.")

    # except Exception as e:
    #     print(f"Error when running ped_analysis: {e}")
    #     exit(1)