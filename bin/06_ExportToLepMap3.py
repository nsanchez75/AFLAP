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

    list_of_Gs = get_LA_info()
    for G_info in list_of_Gs:
        G, LO, UP, P0, SEX = G_info
        print(f"\t Working on parent {G}...")

        forlepmap_file = f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv"
        if os.path.exists(forlepmap_file) and os.path.getsize(forlepmap_file):
            print(f"LepMap3 info for parent {G} has been detected. Skipping.")
            continue

        # determine male and female parent
        sex_dict = {'male': None, 'female': None}
        with open("AFLAP_tmp/Crosses.txt", 'r') as fcrosses:
            for cross in fcrosses:
                cross = cross.strip().split()

                if G == cross[2]:   sex_dict['male'], sex_dict['female'] = G, P0
                elif G == cross[3]: sex_dict['male'], sex_dict['female'] = P0, G
                else: continue
                break

        data = [["CHR", "POS", f"{sex_dict['male']}x{sex_dict['female']}", f"{sex_dict['male']}x{sex_dict['female']}"],
                ["CHR", "POS", sex_dict['male']                          , sex_dict['female']                        ],
                ["CHR", "POS", '0'                                       , '0'                                       ],
                ["CHR", "POS", '0'                                       , '0'                                       ],
                ["CHR", "POS", '1'                                       , '2'                                       ],
                ["CHR", "POS", '0'                                       , '0'                                       ]]
        lepmap_df = pd.DataFrame(data)

        # put all F1 progeny of parent into header
        # TODO: refactor to use df
        with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as fprog1:
            p1_set = set()
            fprog1.readline()
            for p1 in fprog1:
                p1 = p1.strip().split()

                if G in (p1[3], p1[4]) and p1[0] not in p1_set:
                    added_data = [f"{sex_dict['male']}x{sex_dict['female']}", p1[0],
                                    sex_dict['male'], sex_dict['female'], '0', '0']
                    lepmap_df.insert(len(lepmap_df.columns), len(lepmap_df.columns), added_data)
                p1_set.add(p1[0])

        # add rows from filtered tsv file to lepmap df
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
                else: continue
                break

        ftsv = ftsv.replace([0, 1, 2], ['1 0 0 0 0 0 0 0 0 0', '0 1 0 0 0 0 0 0 0 0', '0 0 0 0 1 0 0 0 0 0'], regex=True)

        print(ftsv)
        print(lepmap_df)

        ftsv = ftsv.set_axis(list(lepmap_df.columns), axis=1)
        # add filtered tsv data to lepmap df
        lepmap_df = pd.concat([lepmap_df, ftsv], ignore_index=True)

        # export lepmap df to tsv
        lepmap_df.to_csv(forlepmap_file, sep='\t', header=False, index=False)
        if not os.path.exists(forlepmap_file):
            exit(f"An error occurred: tsv file for {G} has not been created.")
        print(f"\tCompleted making a LepMap3 tsv file for {G}.")
