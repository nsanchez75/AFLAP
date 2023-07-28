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
    for f_type in ["F1", "F2"]:
        if not len(os.listdir(f"AFLAP_tmp/01/{f_type}Count")):
            print(f"No {f_type} progeny found. Skipping.")
            continue

        for G_info in list_of_Gs:
            G, LO, UP, P0, SEX = G_info
            print(f"Working on {f_type} progeny from parent {G}...")

            forlepmap_file = f"AFLAP_Results/{G}_{f_type}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv"
            if os.path.exists(forlepmap_file) and os.path.getsize(forlepmap_file):
                print(f"LepMap3 info for parent {G} has been detected. Skipping.")
                continue

            # determine male and female parent
            sex_dict = {'male': None, 'female': None}
            if SEX == "male":     sex_dict['male'], sex_dict['female'] = G, P0
            elif SEX == "female": sex_dict['male'], sex_dict['female'] = P0, G

            data = [["CHR", "POS", f"{sex_dict['male']}x{sex_dict['female']}", f"{sex_dict['male']}x{sex_dict['female']}"],
                    ["CHR", "POS", sex_dict['male']                          , sex_dict['female']                        ],
                    ["CHR", "POS", '0'                                       , '0'                                       ],
                    ["CHR", "POS", '0'                                       , '0'                                       ],
                    ["CHR", "POS", '1'                                       , '2'                                       ],
                    ["CHR", "POS", '0'                                       , '0'                                       ]]
            lepmap_df = pd.DataFrame(data)

            # put all progeny of parent into header
            prog_df = pd.read_csv(f"AFLAP_tmp/Pedigree_{f_type}.txt", sep='\t')

            prog_df = prog_df[(prog_df["MP"].astype(str) == G) | (prog_df["FP"].astype(str) == G)]
            for prog in prog_df["Individual"].unique():
                added_data = pd.Series([f"{sex_dict['male']}x{sex_dict['female']}", prog,
                            sex_dict['male'], sex_dict['female'], '0', '0'])
                lepmap_df = pd.concat([lepmap_df, added_data], axis=1)

            # add rows from filtered tsv file to lepmap df
            filtered_tsv = f"AFLAP_tmp/05/{G}_{f_type}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv"
            if (not os.path.exists(filtered_tsv)):
                exit(f"An error occurred: Filtered {filtered_tsv} file not found. Rerun 05_ObtainSegStats.py.")
            ftsv = pd.read_csv(filtered_tsv, sep='\t')
            ftsv.insert(0, "MarkerLoc", ftsv["MarkerID"].astype(str) + '_' + ftsv["MarkerLength"].astype(str))
            ftsv = ftsv.drop(["MarkerID", "MarkerLength"], axis=1)

            if SEX == "male":
                ftsv.insert(2, "Male Parent", 1)
                ftsv.insert(3, "Female Parent", 0)
            elif SEX == "female":
                ftsv.insert(2, "Male Parent", 0)
                ftsv.insert(3, "Female Parent", 1)

            ftsv = ftsv.replace([0, 1, 2], ['1 0 0 0 0 0 0 0 0 0', '0 1 0 0 0 0 0 0 0 0', '0 0 0 0 1 0 0 0 0 0'], regex=True)
            ftsv = ftsv.set_axis(list(lepmap_df.columns), axis=1)

            # create lepmap tsv data
            lepmap_df = pd.concat([lepmap_df, ftsv], ignore_index=True)
            lepmap_df.to_csv(forlepmap_file, sep='\t', header=False, index=False)

            if not os.path.exists(forlepmap_file):
                exit(f"An error occurred: {f_type} tsv file for {G} has not been created.")
            print(f"\tCompleted making an {f_type} LepMap3 tsv file for {G}.")
