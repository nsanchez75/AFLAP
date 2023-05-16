import argparse
import os
import subprocess
import sys

from bin.ped_analysis import pedigree_analysis
from bin.helper_funcs.yon import y_or_n

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='AFLAP', description="A script to run all stages of AFLAP.")
    parser.add_argument('-P', '--Pedigree', type=str, required=True, help="Pedigree file (required). See AFLAP README for more information.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for JELLYFISH counting (optional). Default [4].')
    parser.add_argument('-r', '--remove', help='Individual to remove. All other options will be ignored.')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--SDL', type=float, default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--SDU', type=float, default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    parser.add_argument('-k', '--kinship', action='store_true', help='Run kinship estimation.')
    parser.add_argument('-x', '--LowCov', action='store_true', help='Run with low coverage parameters.')
    parser.add_argument('-U', '--Max', type=int, help='Maximum number of markers to output in the genotype tables output under ./AFLAP_Results/')
    args = parser.parse_args()

    # check for dependencies
    print("Checking for dependencies used in AFLAP...")
    for module in ["jellyfish", "ABYSS", "lepmap3"]:
        try:
            subprocess.check_output(args=f"ls $CONDA_PREFIX/bin | grep {module}", shell=True)
        except:
            print(f"Error: {module} not detected.")
            sys.exit(1)
    print("All dependencies found.")

    # make directories if necessary
    os.makedirs("AFLAP_tmp", exist_ok=True)
    os.makedirs("AFLAP_tmp/01", exist_ok=True)

    # perform Pedigree file analysis
    print("Performing pedigree file analysis...")
    if os.path.exists("AFLAP_tmp/PedigreeInfo.txt"):
        with open("AFLAP_tmp/PedigreeInfo.txt") as fped:
            src = fped.readline().strip().split()[1]
            if args.Pedigree == src:
                print(f"\t{args.Pedigree} had already been analyzed. Skipping pedigree analysis.")
            else:
                overwrite_check = y_or_n(f"\t{src}, not {args.Pedigree}, had previously been analyzed. Would you like to overwrite the analysis of {src}? (y/n)")
                if overwrite_check:
                    print(f"\tOverwriting {src}...")
                    pedigree_analysis(args.Pedigree)
                else:
                    print(f"\tNot overwriting {src}. Use it as the pedigree input instead.")
                    sys.exit(0)
    else:
        print(f"\tAnalyzing {args.Pedigree}...")
        pedigree_analysis(args.Pedigree)

    # run programs (#TODO?: allow user to specify what programs to run)
    if args.remove:
        # TODO: figure out what args.remove does
        print("Remove argument passed")
        pass

    print("continue to debugging 01?")
    exit(0)

    DIR = os.path.dirname(os.path.abspath(__file__))

    # TODO: fix all error printing
    # 01_JELLYFISH.py
    try:
        os.system(f"python3 {DIR}/bin/01_JELLYFISH.py -t {args.threads} -m {args.kmer}")
    except SystemExit:
        print("Error found in 01_JELLYFISH.py.")
        exit(1)
    
    # 02_ExtractSingleCopyMers.py
    try:
        os.system(f"python3 {DIR}/bin/02_ExtractSingleCopyMers.py -m {args.kmer}")
    except SystemExit:
        print("Error found in 02_ExtractSingleCopyMers.py.")
    
    # 03_ObtainMarkers.py
    try:
        os.system(f"python3 {DIR}/bin/03_ObtainMarkers.py -m {args.kmer}")
    except SystemExit:
        print("Error found in 03_ObtainMarkers.py.")
    
    # 04_Genotyping.py
    try:
        os.system(f"python3 {DIR}/bin/04_Genotyping.py -m {args.kmer} -L {args.LOD}")
    except SystemExit:
        print("Error found in 04_Genotyping.py.")
    
    # 05_ObtainSegStats.py
    try:
        os.system(f"python3 {DIR}/bin/05_ObtainSegStats.py -m {args.kmer} -L {args.LOD}")
    except:
        print("Error found in 05_ObtainSegStats.py.")
    
    if (args.kinship):
        # TODO: implement 05b
        pass
    elif (args.Max is not None):
        # TODO: implement 05c
        pass

    try:
        os.system(f"python3 {DIR}/bin/06_ExportToLepMap3.py -m {args.kmer}")
    except:
        print("Error found in 06_ExportToLepMap3.py.")

    try:
        os.system(f"python3 {DIR}/bin/07_LepMap3.py -m {args.kmer} -t {args.threads} -L {args.LOD}")
    except:
        print("Error found in 07_LepMap3.py")