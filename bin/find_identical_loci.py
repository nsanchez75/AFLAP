import pandas as pd

# TODO: create 2 separate dataframes by parent and compare them to find equivalent locus sequences

def find_identical_loci(mp_seqs:pd.DataFrame, fp_seqs:pd.DataFrame, mp:str, fp:str)->None:
    # TODO: refactor if we can have more than one parent of same sex
    mp_seqs = mp_seqs.rename(columns={"Sequence": "Male Sequence"})
    fp_seqs = fp_seqs.rename(columns={"Sequence": "Female Sequence"})
    comb_seqs = pd.merge(mp_seqs, fp_seqs, on="Locus Sequence", how='inner')
    comb_seqs.to_csv(f"AFLAP_tmp/03/SimGroups/{mp}_and_{fp}_identical_loci.txt", sep='\t')

# TESTING PURPOSES
# if __name__ == "__main__":
