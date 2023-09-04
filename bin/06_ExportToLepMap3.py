import argparse
import os
import pandas as pd

from get_LA_info import get_LA_info

#################################################
#       A shell script to export the genotype table to LepMap3.
#################################################

def create_f1_forlepmap(kmer:int)->None:
    for G_info in get_LA_info():
        G, LO, UP, P0, SEX = G_info
        print(f"Working on F1 progeny from parent {G}...")

        forlepmap_file = f"AFLAP_Results/{G}_F1_m{kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv"

        # determine male and female parent
        if SEX == "male":     male, female = G, P0
        elif SEX == "female": male, female = P0, G

        # get info from filtered genotype table
        filtered_tsv = f"AFLAP_tmp/05/{G}_F1_m{kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv"
        if (not os.path.exists(filtered_tsv)):
            exit(f"An error occurred: Filtered {filtered_tsv} file not found. Rerun 05_ObtainSegStats.py.")
        ftsv_df = pd.read_csv(filtered_tsv, sep='\t')
        ftsv_df.insert(0, "MarkerLoc", ftsv_df["MarkerID"].astype(str))
        ftsv_df = ftsv_df.drop(["MarkerID"], axis=1)

        # add columns for male and female parents
        if SEX == "male":
            ftsv_df.insert(2, "Male Parent", 1)
            ftsv_df.insert(3, "Female Parent", 0)
        elif SEX == "female":
            ftsv_df.insert(2, "Male Parent", 0)
            ftsv_df.insert(3, "Female Parent", 1)

        ftsv_df = ftsv_df.replace([0, 1, 2], ['1 0 0 0 0 0 0 0 0 0', '0 1 0 0 0 0 0 0 0 0', '0 0 0 0 1 0 0 0 0 0'], regex=False)

        # create lepmap header
        data = [["CHR", "POS", f"{male}x{female}", f"{male}x{female}"],
                ["CHR", "POS", male              , female            ],
                ["CHR", "POS", '0'               , '0'               ],
                ["CHR", "POS", '0'               , '0'               ],
                ["CHR", "POS", '1'               , '2'               ],
                ["CHR", "POS", '0'               , '0'               ]]
        lepmap_df = pd.DataFrame(data)

        # put all progeny of parent into header
        prog_df = pd.read_csv(f"AFLAP_tmp/Pedigree_F1.txt", sep='\t')
        prog_df = prog_df[(prog_df["MP"].astype(str).isin(male.split('_'))) | (prog_df["FP"].astype(str).isin(female.split('_')))]
        included_progs = set(ftsv_df.columns[4:])
        for prog in prog_df["Individual"].unique():
            if prog in included_progs:
                added_data = pd.Series([f"{male}x{female}", prog, male, female, '0', '0'])
                lepmap_df = pd.concat([lepmap_df, added_data], axis=1)

        # create lepmap data file
        ftsv_df = ftsv_df.set_axis(list(lepmap_df.columns), axis=1)
        lepmap_df = pd.concat([lepmap_df, ftsv_df], ignore_index=True)
        lepmap_df.to_csv(forlepmap_file, sep='\t', header=False, index=False)

        if not os.path.exists(forlepmap_file):
            exit(f"An error occurred: F1 tsv file for {G} has not been created.")
        print(f"\tCompleted making an F1 LepMap3 tsv file for {G}.")

def create_f2_forlepmap(kmer:int)->None:
    male = list()
    female = list()
    for G, _, _, _, SEX in get_LA_info():
        if SEX == "male": male.append(G)
        else:             female.append(G)
    male = '_'.join(male)
    female = '_'.join(female)

    forlepmap_file = f"AFLAP_Results/{male}x{female}_F2_m{kmer}.ForLepMap3.tsv"

    # get info from filtered genotype table
    filtered_tsv = f"AFLAP_tmp/05/{male}x{female}_F2_m{kmer}.Genotypes.MarkerID.Filtered.tsv"
    if (not os.path.exists(filtered_tsv)):
        exit(f"An error occurred: Filtered {filtered_tsv} file not found. Rerun 05_ObtainSegStats.py.")
    ftsv_df = pd.read_csv(filtered_tsv, sep='\t')
    ftsv_df = ftsv_df.iloc[:, [1, 0] + list(range(2, ftsv_df.shape[1]))]

    # add columns for male and female parents (F0 and F1)
    ftsv_df.insert(2, "Male F0 Parent", 'AA')
    ftsv_df.insert(3, "Female F0 Parent", 'BB')
    ftsv_df.insert(4, "Male F1 Parent", 'AB')
    ftsv_df.insert(5, "Female F1 Parent", 'AB')

    # replace genotype table values to their chromosome data equivalents
    ftsv_df = ftsv_df.replace(['XX', 'AA', 'BB', 'AB'], ['0.33 0.33 0 0 0.33 0 0 0 0 0', '1 0 0 0 0 0 0 0 0 0', '0 0 0 0 1 0 0 0 0 0', '0 1 0 0 0 0 0 0 0 0'], regex=False)

    # create lepmap header
    data = [["CHR", "POS", f"{male}x{female}" , f"{male}x{female}", f"{male}x{female}", f"{male}x{female}"],
            ["CHR", "POS", male               , female            , "DUM1"            , "DUM2"            ],
            ["CHR", "POS", '0'                , '0'               , male              , male              ],
            ["CHR", "POS", '0'                , '0'               , female            , female            ],
            ["CHR", "POS", '1'                , '2'               , '1'               , '2'               ],
            ["CHR", "POS", '0'                , '0'               , '0'               , '0'               ]]
    lepmap_df = pd.DataFrame(data)

    # put all progeny into header
    prog_df = pd.read_csv(f"AFLAP_tmp/Pedigree_F2.txt", sep='\t')
    prog_df = prog_df[(prog_df["MP"].astype(str).isin(male.split('_'))) | (prog_df["FP"].astype(str).isin(female.split('_')))]
    included_progs = set(ftsv_df.columns[6:])
    for prog in prog_df["Individual"].unique():
        if prog in included_progs:
            added_data = pd.Series([f"{male}x{female}", prog, "DUM1", "DUM2", '0', '0'])
            lepmap_df = pd.concat([lepmap_df, added_data], axis=1)

    # create lepmap data file
    ftsv_df = ftsv_df.set_axis(list(lepmap_df.columns), axis=1)
    lepmap_df = pd.concat([lepmap_df, ftsv_df], ignore_index=True)
    lepmap_df.to_csv(forlepmap_file, sep='\t', header=False, index=False)

    if not os.path.exists(forlepmap_file):
        exit(f"An error occurred: F2 tsv file for {male}x{female} has not been created.")
    print(f"\tCompleted making an F2 LepMap3 tsv file for {male}x{female}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ExportToLepMap3', description="A script to export the genotype table to LepMap3.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()

    for f_type in ["F1", "F2"]:
        if not len(os.listdir(f"AFLAP_tmp/01/{f_type}Count")):
            print(f"No {f_type} progeny found. Skipping.")
            continue

        if f_type == "F1": create_f1_forlepmap(args.kmer)
        else:              create_f2_forlepmap(args.kmer)
