import argparse
import os
import subprocess

import bin.check_newfiles as check_nf
import bin.pedigree_check as pedigree_check


def main():
    parser = argparse.ArgumentParser(prog='AFLAP', description='A script to run all stages of AFLAP.')
    parser.add_argument('-P', '--Pedigree', required=True, help='Pedigree file (required). See AFLAP README for more information.')
    parser.add_argument('-m', '--kmer', default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-t', '--threads', default=4, help='Threads for JELLYFISH counting (optional). Default [4].')
    parser.add_argument('-r', '--remove', help='Individual to remove. All other options will be ignored.')
    parser.add_argument('-L', '--LOD', default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--Sdl', default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--Sdu', default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    parser.add_argument('-k', '--kinship', action='store_true', help='Run kinship estimation.')
    parser.add_argument('-x', '--LowCov', action='store_true', help='Run with low coverage parameters.')
    parser.add_argument('-U', '--Max', help='Maximum number of markers to output in the genotype tables output under ./AFLAP_Results/')
    args = parser.parse_args()

    # argument check
    if args.threads is None:
        print("Threads not specified, will proceed with default [4]")
    if args.kmer is None:
        print("mer size not specified, will proceed with default [31]")


    # jellyfish check
    try:
        jellyfish_version = subprocess.check_output("jellyfish --version", shell=True)
        jellyfish_version = jellyfish_version.decode().strip()
        print(f"{jellyfish_version} detected\n")
    except OSError:
        print("jellyfish not detected, please modify your PATH. Terminating.")
        exit(1)


    # make a temporary directory
    os.makedirs("AFLAP_tmp", exist_ok=True)
    # strip '#' from pedigree file and make new temporary file for it
    with open(args.Pedigree, 'r') as infile, open("AFLAP_tmp/Pedigree.txt", 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                outfile.write(line)
    
    args.Pedigree = "AFLAP_tmp/Pedigree.txt"

    # separate individuals in pedigree file
    pedigree_check.separate_ped("AFLAP_tmp/Pedigree.txt")
    # check if F1 and F2 files are formatted correctly
    if os.path.exists("AFLAP_tmp/Pedigree_F1.txt"):
        pedigree_check.check_format("AFLAP_tmp/Pedigree_F1.txt")
    elif os.path.exists("AFLAP_tmp/Pedigree_F2.txt"):
        pedigree_check.check_format("AFLAP_tmp/Pedigree_F2.txt")
    
    # make 01 directory
    os.makedirs("AFLAP_tmp/01", exist_ok=True)
    # check the Pedigree file so that parents are specified
    pedigree_check.confirm_parents()
    # check pedigree structure (Currently only tested on F1 and F2 independently. Has not been tested on both simultaneously.)
    pedigree_check.check_pedigree_structure()   # TODO: get rid of F2 restrictions from this function eventually


    if args.remove is None:     # TODO: understand what args.remove is
        DIR = os.path.dirname(os.path.abspath(__file__))        
        
        os.system(f'python3 {DIR}/bin/01_JELLYFISH.py -P {args.Pedigree} -t {args.threads} -m {args.kmer}')
        check_nf.check_01_jellyfish()
        
        os.system(f'python3 {DIR}/bin/02_ExtractSingleCopyMers.py -P {args.Pedigree} -m {args.kmer}')
        check_nf.check_02_escm()

        os.system(f'python3 {DIR}/bin/03_ObtainMarkers.py -P {args.Pedigree} -m {args.kmer}')

        if args.LowCov:
            os.system(f'{DIR}/bin/04_Genotyping.sh -P {args.Pedigree} -m {args.kmer} -L')
            
            print("continue to 05?")
            exit(0)

            os.system(f'{DIR}/bin/05_ObtainSegStats.sh -d {args.Sdl} -D {args.Sdu} -P {args.Pedigree} -m {args.kmer} -L')
        else:
            os.system(f'{DIR}/bin/04_Genotyping.sh -P {args.Pedigree} -m {args.kmer}')
            os.system(f'{DIR}/bin/05_ObtainSegStats.sh -d {args.Sdl} -D {args.Sdu} -P {args.Pedigree} -m {args.kmer}')
        if args.kinship:
            os.system(f'{DIR}/bin/06_KinshipEstimation.sh -P {args.Pedigree}')
        if args.Max:
            os.system(f'{DIR}/bin/07_OutputMarkers.sh -P {args.Pedigree} -m {args.kmer} -U {args.Max}')

if __name__ == '__main__':
    main()
