import pandas as pd

from get_LA_info import get_LA_info

def find_identical_loci(seq_groups:pd.DataFrame)->None:
    print("Finding identical sequence loci...")

    # identify parents categorized by sex
    crosses_df = pd.read_csv("AFLAP_tmp/Crosses.txt", sep='\t', usecols=[2, 3], names=["MP", "FP"])
    MP = set(crosses_df["MP"].unique())
    FP = set(crosses_df["FP"].unique())

    # create a dataframe of same loci sequences
    seqs_df = pd.DataFrame(columns=["Sequence"])
    for useq in seq_groups["Locus Sequence"].unique():
        useq_df = seq_groups[seq_groups["Locus Sequence"] == useq]
        if len(useq_df.index) < 2: continue
        print(useq_df)
        parents = useq_df["Parent"].unique().tolist()
        # detect if both parents identified to locus sequence
        mpcheck = False
        fpcheck = False
        for p in parents:
            if p in MP: mpcheck = True
            if p in FP: fpcheck = True
        # if identified add sequences associated to locus sequence to tsv file dataframe
        if mpcheck and fpcheck:
            for seq in useq_df["Sequence"].unique():
                if seq not in set(seqs_df["Sequence"].unique()):
                    seqs_df.loc[len(seqs_df.index)] = [seq]

    # produce a tsv file of same-loci sequences
    seqs_df.to_csv("AFLAP_tmp/03/SameLociSeqs.tsv", sep='\t', index=False, header=False)

# TESTING PURPOSES
if __name__ == "__main__":
    seq_groups = pd.DataFrame(columns=["Sequence", "Locus Sequence", "Parent"])
    kmer = 31
    ak = 2 * kmer - 1

    print("Running test on find_identical_loci.py")

    for G_info in get_LA_info():
        G, LO, UP, P0 = G_info

        print(f"Analyzing sequences from parent {G}...")
        with open(f"AFLAP_tmp/03/{G}_m{kmer}_L{LO}_U{UP}_abyss.fa", 'r') as fab:
            for i in range(500):
                id = fab.readline().strip()
                seq = fab.readline().strip()
                if not id or not seq: break

                id = id.replace('>', '').split()
                if int(id[1]) >= ak:
                    subseq = seq[9:(9 + kmer)]
                    rc_subseq = subseq[::-1].translate(subseq.maketrans("ATCG", "TAGC"))
                    if subseq > rc_subseq: subseq = rc_subseq

                    seq_groups.loc[len(seq_groups.index)] = [subseq, seq[0:(kmer - 1)] + seq[(len(seq) - kmer + 1):], G]

        print(f"Finished analyzing sequences from parent {G}.")

    print("Sequence groups:")
    print(seq_groups)

    print("Finding identical loci...")
    find_identical_loci(seq_groups)
    print("Test complete")