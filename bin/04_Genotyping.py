import argparse
import os


def main():
    parser = argparse.ArgumentParser(prog='Genotyping', description="A script to genotype progeny")
    parser.add_argument('-P', '--Pedigree', required=True, help='Pedigree file (required). See AFLAP README for more information.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()


    # make directories
    os.makedirs(["AFLAP_tmp/04/Count", "AFLAP_tmp/04/Call"], exist_ok=True)


if __name__ == "__main__":
    main()