import pandas as pd

def find_identifier(parents_set:set, useq_df:pd.DataFrame)->str:
    identifier = useq_df["Identifier"][useq_df["Parent"].isin(parents_set)].tolist()
    if len(identifier) != 1: raise ValueError("Invalid number of parents of the same sex for one individual.")
    return identifier[0]

def find_identical_loci(seq_groups:pd.DataFrame)->None:
    print("Finding identical sequence loci...")

    # identify parents categorized by sex
    crosses_df = pd.read_csv("AFLAP_tmp/Crosses.txt", sep='\t', usecols=[2, 3], names=["MP", "FP"])
    MP = set(crosses_df["MP"].unique())
    FP = set(crosses_df["FP"].unique())

    # create a dataframe of same loci sequences
    seqs_df = pd.DataFrame(columns=["Sequence", "Male Identifier", "Female Identifier"])
    for useq in seq_groups["Sequence"].unique():
        useq_df = seq_groups[seq_groups["Sequence"] == useq]
        male_identifier = find_identifier(MP, useq_df)
        female_identifier = find_identifier(FP, useq_df)
        seqs_df.loc[len(seqs_df.index)] = [useq, male_identifier, female_identifier]

    # produce a tsv file of same-loci sequences
    seqs_df.to_csv("AFLAP_tmp/03/SameLociSeqs.tsv", sep='\t', index=False)