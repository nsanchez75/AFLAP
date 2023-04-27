import argparse
import os
import subprocess
import sys

import bin.ped_analysis as pa

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


    # 1. check for modules
    try:
        jellyfish_version = subprocess.check_output("jellyfish --version", shell=True).decode().strip()
        print(f"{jellyfish_version} detected.")
    except OSError:
        print("Error in AFLAP.py: Jellyfish not detected. Please modify your PATH.")
        sys.exit(1)
    try:
        abyss_version = subprocess.check_output("ABYSS --version", shell=True).decode().strip()
        print(f"{abyss_version} detected.")
    except OSError:
        print("Error in AFLAP.py: ABySS not detected. Please modify your PATH.")
        sys.exit(1)
    # TODO: add R/4.0.1 as well?


    # 2. perform Pedigree file analysis
    # if os.path.exists("AFLAP_tmp/Pedigree.txt"):
    #     # TODO: check if file is same as argument input
    #     pass
    # else:
    os.makedirs("AFLAP_tmp", exist_ok=True)
    os.makedirs("AFLAP_tmp/01", exist_ok=True)

    pa.pedigree_analysis(args.Pedigree)


    # 3. run programs (#TODO?: allow user to specify what programs to run)
    if args.remove:
        # TODO: figure out what args.remove does
        print("Remove argument passed")
        pass

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