import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(prog='Genotyping', description="A script to genotype progeny")
    parser.add_argument('-P', '--Pedigree', required=True, help='Pedigree file (required). See AFLAP README for more information.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-L', '--LOD', default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()


    # make directories
    os.makedirs(["AFLAP_tmp/04/Count", "AFLAP_tmp/04/Call"], exist_ok=True)


    # check for markers
    with open("AFLAP_tmp/01.LA.txt", 'r') as fla:
        for g in fla:
            g = g.strip()
            print(f"Calling GT for {g} derived markers...")

            # identify correct file based on intermediate results
            with open(f"AFLAP_tmp/03/{g}_CrossedTo.txt", 'r') as fcrossed, open("AFLAP_tmp/02/Boundaries.txt", 'r') as fbounds:
                P0 = fcrossed.read().strip()

                for bounds in fbounds:
                    bounds = bounds.strip().split()

                    if bounds[0] == g:
                        LO = bounds[1]
                        HI = bounds[2]

            if os.path.exists(f"AFLAP_tmp/03/ParentalMarkers/{g}_m{args.kmer}_MARKERS_L{LO}_U{HI}_{P0}.fa"):
                with open(f"AFLAP_tmp/03/ParentalMarkers/{g}_m{args.kmer}_MARKERS_L{LO}_U{HI}_{P0}.fa", 'r') as fmark:
                    m_count = 0
                    for line in fmark:
                        if line.startswith('>'):
                            m_count += 1

                    print(f"{m_count} markers identified in AFLAP_tmp/03/ParentalMarkers/{g}_m{args.kmer}_MARKERS_L{LO}_U{HI}_{P0}.fa. These will be surveyed against progeny.")
            else:
                print("Can't automatically locate marker file. Please try rerunning AFLAP.py.")
                exit(1)     # FIXME: doesn't exit in .sh file check on this
            
            # get lines in pedigree file whose parent is g
            h_list = []
            if os.path.exists("AFLAP_tmp/Pedigree_F1.txt"):
                with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as ff1:
                    for line in ff1:
                        line = line.strip().split()

                        if g in {line[3], line[4]}:
                            h_list += [line]
            if os.path.exists("AFLAP_tmp/Pedigree_F2.txt"):
                with open("AFLAP_tmp/Pedigree_F2.txt", 'r') as ff2:
                    for line in ff2:
                        line = line.strip().split()

                        if g in {line[3], line[4]}:
                            h_list += [line]
            
            if not len(h_list):
                print(f"Error: No progeny found whose parent is {g}. Terminating.")
                exit(1)
            
            for h in h_list:
                if not os.path.exists(f"AFLAP_tmp/01/ProgCounts/{h}.jf{args.kmer}"):
                    print("Error: No has for {h} detected. Please rerun 01_JELLYFISH.py.")
                    exit(1)

                if os.path.exists(f"AFLAP_tmp/04/Count/{h}_{g}_m{args.kmer}_L{LO}_U{HI}_{P0}.txt") and \
                   os.path.exists(f"AFLAP_tmp/04/Call/{h}_{g}_m{args.kmer}_L{LO}_U{HI}_{P0}.txt"):
                    print(f"Genotype of {g} markers for {h} previously calculated, will use these.")
                else:
                    jellyfish_cmd = f"jellyfish query -s AFLAP_tmp/03/ParentalMarkers/{g}_m{args.kmer}_MARKERS_L{LO}_U{HI}_{P0}.fa AFLAP_tmp/01/ProgCounts/{h}.jf{args.kmer} > AFLAP_tmp/04/Count/{h}_{g}_m{args.kmer}_L{LO}_U{HI}_{P0}.txt"
                    jellyfish_output = subprocess.run(jellyfish_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout.split('\n')
                    with open(f"AFLAP_tmp/04/Count/{h}_{g}_m{args.kmer}_L{LO}_U{HI}_{P0}.txt", 'r') as fin, open(f"AFLAP_tmp/04/Call/{h}_{g}_m{args.kmer}_L{LO}_U{HI}_{P0}.txt", 'w') as fout:
                        for line in fin:
                            line = line.strip().split()

                            if int(line[1]) >= args.LOD:
                                fout.write("1\n")
                            else:
                                fout.write("0\n")
            
            print(f"GT calling for {g} derived markers complete!")
            
            


if __name__ == "__main__":
    main()