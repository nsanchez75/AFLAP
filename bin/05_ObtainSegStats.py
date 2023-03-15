import argparse
import os
import sys

def main()->None:
    parser = argparse.ArgumentParser(prog='ObtainSegStats', description='A script to plot marker distributions in progeny.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()


    # 1. create directories
    os.makedirs("AFLAP_tmp/05/FilteredCall", exist_ok=True)
    os.makedirs("AFLAP_tmp/05/SegregationInfo", exist_ok=True)


    # 2. initialize variables
    ak = 2 * int(args.kmer) - 1


    # 3. analyze all parents whom we can create a genetic map for
    print("Performing segment statistics analysis...")
    with open("AFLAP_tmp/01/LA.txt", 'r') as fla:
        for p in fla:
            # extract info about analyzed parent
            p = p.strip().split()
            G = p[0]
            LO = p[1]
            UP = p[2]

            print(f"\tObtaining segment statistics for {G}...")

            with open(f"AFLAP_tmp/03/{G}_CrossedTo.txt", 'r') as fct:
                p0 = []
                for op in fct:
                    p0.append(op.strip())
                p0 = '_'.join(p0)

            # check for genotype table
            if not os.path.exists(f"AFLAP_tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.MarkerID.tsv"):
                print(f"Error in 05_ObtainSegStats.py: Genotype table for {G} not found. Rerun 04_Genotyping.py.")
                sys.exit(1)

            # count progeny
            with open(f"AFLAP/tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.Marker.tsv", 'r') as ftsv:
                num_prog = len(ftsv.readline()) - 2

            print(f"\t\t{num_prog} Genotype calls for {G} detected. Summarizing...")

            # compare contig lengths that subsequences came from TODO: fact check whether it comes from ABySS or jellyfish
            with open(f"AFLAP/tmp/04/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}.Genotypes.Marker.tsv", 'r') as ftsv:
                # initialize dictionaries for MarkerEqual and MarkerOver
                meq = dict()
                mov = dict()

                for tsv_line in ftsv:
                    # extract original contig length
                    tsv_line = tsv_line.strip().split()
                    clen = tsv_line[1].split('_')[1]

                    # initialize number of calls variable
                    num_calls = 0

                    # count number of calls for sequence in progeny
                    for i in range(2, len(tsv_line)):
                        num_calls += int(tsv_line[i])

                    # find proportion of calls to progeny
                    prop_calls = num_calls / num_prog

                    # add values to respective dictionaries
                    if clen == ak:
                        if prop_calls not in meq:
                            meq[prop_calls] = 1
                        else:
                            meq[prop_calls] += 1
                    else:
                        if prop_calls not in mov:
                            mov[prop_calls] = 1
                        else:
                            meq[prop_calls] += 1

                # write to MarkerEqual and MarkerOver
                with open(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerEqual{ak}.hist", 'a') as fme:
                    for m in meq:
                        fme.write(f"{m} {meq[m]}\n")
                with open(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{p0}_MarkerOver{ak}.hist", 'a') as fmo:
                    for m in mov:
                        fmo.write(f"{m} {mov[m]}\n")


if __name__ == "__main__":
    main()