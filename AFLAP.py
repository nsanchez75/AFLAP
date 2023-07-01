import argparse
import os
import subprocess

from bin.ped_analysis import pedigree_analysis
from bin.marker_reduction import marker_reduction
from bin.update_individual import update_individual

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='AFLAP', description="A script to run all stages of AFLAP.")
    parser.add_argument('-P', '--Pedigree', type=str, required=True, help="Pedigree file (required). See AFLAP README for more information.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for JELLYFISH counting (optional). Default [4].')
    parser.add_argument('-r', '--remove', type=str, help='Individual to remove. All other options will be ignored.')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--SDL', type=float, default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--SDU', type=float, default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    parser.add_argument('-x', '--LowCov', type=int, default=2, help='Run with low coverage parameters.')
    parser.add_argument('-U', '--Max', type=int, help='Maximum number of markers to output in the genotype tables output under ./AFLAP_Results/')
    args = parser.parse_args()

    # run update function if necessary
    if args.remove:
        print("Performing individual removal...")
        update_individual(args.remove, args.Pedigree)
        print(f"Successfully removed individual {args.remove}. Rerun AFLAP to regenerate the removed files.")
        exit(0)

    # check for dependencies
    print("Checking for dependencies used in AFLAP...")
    for module in ["jellyfish", "ABYSS", "lepmap3"]:
        try:
            subprocess.check_output(args=f"ls $CONDA_PREFIX/bin | grep {module}", shell=True)
        except:
            print(f"Error: {module} not detected.")
            exit(1)
    print("All dependencies found.")

    # make directories if necessary
    os.makedirs("AFLAP_tmp", exist_ok=True)
    os.makedirs("AFLAP_Results", exist_ok=True)

    # perform Pedigree file analysis
    print("Performing pedigree file analysis...")
    if os.path.exists("AFLAP_tmp/PedigreeInfo.txt"):
        with open("AFLAP_tmp/PedigreeInfo.txt") as fped:
            src = fped.readline().strip().split()[1]
            print(f"\tWARNING: You are overwriting the analysis of {src}.")
            print(f"\tOverwriting {src}...")
    print(f"\tAnalyzing {args.Pedigree}...")
    pedigree_analysis(args.Pedigree)

    # # run programs (#TODO?: allow user to specify what programs to run)
    # if args.remove:
    #     # TODO: figure out what args.remove does
    #     print("Remove argument passed")
    #     pass

    DIR = os.path.dirname(os.path.abspath(__file__))

    try:
        # 01_JELLYFISH.py
        os.system(f"python3 {DIR}/bin/01_JELLYFISH.py -t {args.threads} -m {args.kmer}")

        # 02_ExtractSingleCopyMers.py
        os.system(f"python3 {DIR}/bin/02_ExtractSingleCopyMers.py -m {args.kmer}")

        # 03_ObtainMarkers.py
        os.system(f"python3 {DIR}/bin/03_ObtainMarkers.py -m {args.kmer}")

        # 04_Genotyping.py
        os.system(f"python3 {DIR}/bin/04_Genotyping.py -m {args.kmer} -x {args.LowCov}")

        # 05_ObtainSegStats.py
        os.system(f"python3 {DIR}/bin/05_ObtainSegStats.py -m {args.kmer} -L {args.LOD}")

        if (args.Max is not None): marker_reduction(args.kmer, args.Max)

        # 06_ExportToLepMap3.py
        os.system(f"python3 {DIR}/bin/06_ExportToLepMap3.py -m {args.kmer}")

        # 07_LepMap3.py
        os.system(f"python3 {DIR}/bin/07_LepMap3.py -m {args.kmer} -t {args.threads} -L {args.LOD}")

        print("AFLAP complete!")
    except:
        exit(1)
