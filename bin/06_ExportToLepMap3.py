import argparse
import os
import pandas as pd

from get_LA_info import get_LA_info

#################################################
#       A shell script to export the genotype table to LepMap3.
#################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ExportToLepMap3', description="A script to export the genotype table to LepMap3.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()

    # get parents to run LepMap3 on
    list_of_Gs = get_LA_info()
    for G_info in list_of_Gs:
        G, LO, UP, P0 = G_info

        # print init statement
        print(f"\t Working on parent {G}...")

        if os.path.exists(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv") and os.path.getsize(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv"):
            print(f"LepMap3 info for parent {G} has been detected. Skipping.")
        else:
            # determine male and female parent
            sex_dict = {'male': None, 'female': None}
            with open("AFLAP_tmp/Crosses.txt", 'r') as fcrosses:
                for cross in fcrosses:
                    cross = cross.strip().split()

                    # continue if parent not found
                    if G not in (cross[2], cross[3]): continue

                    # identify sex if parent found
                    if G == cross[2]:
                        sex_dict['male']    = G
                        sex_dict['female']  = P0
                    elif G == cross[3]:
                        sex_dict['male']    = P0
                        sex_dict['female']  = G
                    break

            data =  [["CHR", "POS", f"{sex_dict['male']}x{sex_dict['female']}", f"{sex_dict['male']}x{sex_dict['female']}"],
                        ["CHR", "POS", sex_dict['male']                          , sex_dict['female']                        ],
                        ["CHR", "POS", '0'                                       , '0'                                       ],
                        ["CHR", "POS", '0'                                       , '0'                                       ],
                        ["CHR", "POS", '1'                                       , '2'                                       ],
                        ["CHR", "POS", '0'                                       , '0'                                       ]]
            df = pd.DataFrame(data)

            f1_progs_df = pd.read_csv("AFLAP_tmp/Pedigree_F1.txt", sep='\t')
            f1_progs = f1_progs_df["Individual"].unique().tolist()
            p1_set = set()

            for _, row in f1_progs_df.iterrows():
                p1 = row.to_dict()

                if G in (p1["MP"], p1["FP"]) and p1["Individual"] not in p1_set:
                    added_data = [f"{sex_dict['male']}x{sex_dict['female']}", p1["Individual"],
                                    sex_dict['male'], sex_dict['female'], '0', '0']
                    df.insert(len(df.columns), len(df.columns), added_data)

                p1_set.add(p1["Individual"])

            # put all F1 progeny of parent G into data header
            # with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as fprog1:
            #     p1_set = set()
            #     fprog1.readline()   # TODO: delete later after finished refactoring to implement df
            #     for p1 in fprog1:
            #         p1 = p1.strip().split()

            #         if G in (p1[3], p1[4]) and p1[0] not in p1_set:
            #             added_data = [f"{sex_dict['male']}x{sex_dict['female']}", p1[0],
            #                             sex_dict['male'], sex_dict['female'], '0', '0']
            #             df.insert(len(df.columns), len(df.columns), added_data)

            #         p1_set.add(p1[0])

            # add rows from filtered tsv file to df
            if (not os.path.exists(f"AFLAP_tmp/05/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv")):
                exit("An error occurred: Filtered .tsv file not found. Rerun 05_ObtainSegStats.py.")
            ftsv = pd.read_csv(f"AFLAP_tmp/05/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv", sep='\t')
            ftsv.insert(0, "MarkerLoc", ftsv["MarkerID"].astype(str) + '_' + ftsv["MarkerLength"].astype(str))
            ftsv = ftsv.drop(["MarkerID", "MarkerLength"], axis=1)

            with open("AFLAP_tmp/Crosses.txt", 'r') as fcrosses:
                for cross in fcrosses:
                    cross = cross.strip().split()
                    if (cross[2] == G):
                        ftsv.insert(2, "Male Parent", 1)
                        ftsv.insert(3, "Female Parent", 0)
                    elif (cross[3] == G):
                        ftsv.insert(2, "Male Parent", 0)
                        ftsv.insert(3, "Female Parent", 1)
                    else:
                        continue
                    break

            ftsv = ftsv.replace([0, 1], ['1 0 0 0 0 0 0 0 0 0', '0 1 0 0 0 0 0 0 0 0'], regex=True)
            ## set columns of ftsv DataFrame to match df's column names
            ftsv = ftsv.set_axis(list(df.columns), axis=1)
            ## concatenate ftsv under df
            df = pd.concat([df, ftsv], ignore_index=True)

            # export df DataFrame to tsv dedicated to LepMap3
            df.to_csv(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv", sep='\t', header=False, index=False)
            if not os.path.exists(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv"):
                exit(f"An error occurred: tsv file for {G} has not been created.")

            print(f"\tCompleted making a LepMap3 tsv file for {G}.")
