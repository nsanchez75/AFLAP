import argparse
import os
import pandas as pd
import numpy as np

import get_LA_info as gli

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ExportToLepMap3', description="A script to export the genotype table to LepMap3.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    args = parser.parse_args()


    # 1. make directories
    os.makedirs("AFLAP_tmp/06", exist_ok=True)

    # 2. get parents to run LepMap3 on
    if (not os.path.exists("AFLAP_tmp/01/LA.txt")):
        raise FileNotFoundError("AFLAP_tmp/01/LA.txt not found. Rerun pipeline.")
    list_of_Gs = gli.get_LA_info("AFLAP_tmp/01/LA.txt", "AFLAP_tmp/01/Crosses.txt")
    for G_info in list_of_Gs:
        G, LO, UP, P0 = G_info

        # determine male and female parent
        sex_dict = {'male': None, 'female': None}
        with open("AFLAP_tmp/01/Crosses.txt", 'r') as fcrosses:
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

        # put all F1 progeny of parent G into data header
        with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as fprog1:
            p1_set = set()
            for p1 in fprog1:
                p1 = p1.strip().split()

                if G in (p1[3], p1[4]) and p1[0] not in p1_set:
                    added_data = [f"{sex_dict['male']}x{sex_dict['female']}", p1[0],
                                    sex_dict['male'], sex_dict['female'], '0', '0']
                    df.insert(len(df.columns), len(df.columns), added_data)

                p1_set.add(p1[0])

        # add rows from filtered tsv file to df
        if (not os.path.exists(f"AFLAP_tmp/05/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv")):
            raise FileNotFoundError("Filtered .tsv file not found. Rerun 05_ObtainSegStats.py.")
        ftsv = pd.read_csv(f"AFLAP_tmp/05/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.Genotypes.MarkerID.Filtered.tsv", sep='\t')
        ftsv.insert(1, "MarkerLoc", ftsv["MarkerSequence"].astype(str) + '_' + ftsv["MarkerValue"].astype(str))
        ftsv = ftsv.drop(["MarkerSequence", "MarkerValue"], axis=1)
        ftsv.insert(2, "Male Parent", 1)
        ftsv.insert(3, "Female Parent", 0)
        ftsv = ftsv.replace([0, 1], ['1 0 0 0 0 0 0 0 0 0', '0 1 0 0 0 0 0 0 0 0'], regex=True)

        for i in range(0, len(ftsv.index)):
            df.loc[len(df.index) + i] = ftsv.loc[i].to_numpy()

        print(df)

        # TODO:
        #   - iterate through filtered genotype.tsv's progeny [M1 ... M99]
        #   - convert 0 -> 1 0 0 0 0 0 0 0 0 0, 1 -> 0 1 0 0 0 0 0 0 0 0
        #   - place into data as column