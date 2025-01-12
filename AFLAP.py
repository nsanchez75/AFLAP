import argparse
import os
import subprocess

from bin.ped_analysis import pedigree_analysis
from bin.marker_reduction import marker_reduction
from bin.update_individual import update_individual

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='AFLAP', description="A script to run all stages of AFLAP.")
    parser.add_argument('-P', '--Pedigree', type=str, required=True, help="Pedigree file (required). See AFLAP README for more information.")
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size. Default [31].')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for JELLYFISH counting. Default [4].')
    parser.add_argument('-r', '--remove', type=str, help='Individual to remove. All other options will be ignored.')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-d', '--SDL', type=float, default=0.2, help='Lower boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.2].')
    parser.add_argument('-D', '--SDU', type=float, default=0.8, help='Upper boundary for marker cut off. Can be used to filter for segregation distortion. Default [0.8].')
    parser.add_argument('-f', '--fXX', type=float, default=None, help='Limit for how many XX can exist in a row. If surpassed then sequence is not considered for analysis. Default [None].')
    parser.add_argument('-x', '--LowCov', type=int, default=2, help='Run with low coverage parameters.')
    parser.add_argument('-n', '--nLG', type=int, default=10, help='Minimum number of linkage groups needed to continue F2 LepMap3 analysis. Default [10].')
    parser.add_argument('-U', '--Max', type=int, help='Maximum number of markers to output in the genotype tables output under ./AFLAP_Results/')
    args = parser.parse_args()

    # check for dependencies
    print("Checking for dependencies used in AFLAP...")
    for module in ["jellyfish", "ABYSS", "lepmap3"]:
        try:
            subprocess.check_output(args=f"ls $CONDA_PREFIX/bin | grep {module}", shell=True)
        except:
            exit(f"An error occurred: {module} not detected.")
    print("All dependencies found.")

    # run update function if necessary
    if args.remove:
        print("Performing individual removal...")
        update_individual(args.remove, args.Pedigree)
        print(f"Successfully removed individual {args.remove}. Rerun AFLAP to regenerate the removed files.")
        exit(0)

    # make directories if necessary
    os.makedirs("AFLAP_tmp", exist_ok=True)
    os.makedirs("AFLAP_Results", exist_ok=True)

    # perform Pedigree file analysis
    print("Performing pedigree file analysis...")
    if os.path.exists("AFLAP_tmp/PedigreeInfo.txt"):
        with open("AFLAP_tmp/PedigreeInfo.txt") as fped:
            src = fped.readline().strip().split()[1]
            print(f"\tWARNING: You are overwriting the analysis of {src}.")
    print(f"\tAnalyzing {args.Pedigree}...")
    pedigree_analysis(args.Pedigree)
    print("Information from pedigree file extracted.")

    DIR = os.path.dirname(os.path.abspath(__file__))

    try:
        # 01_JELLYFISH.py
        print("\nStep 1/6: Jellyfish Counting\n")
        subprocess.run(f"python3 {DIR}/bin/01_JELLYFISH.py -t {args.threads} -m {args.kmer}",
                       check=True, shell=True)
        # 02_ExtractSingleCopyMers.py
        print("\nStep 2/6: Extracting Single-Copy K-Mers\n")
        subprocess.run(f"python3 {DIR}/bin/02_ExtractSingleCopyMers.py -m {args.kmer}",
                       check=True, shell=True)
        # 03_ObtainMarkers.py
        print("\nStep 3/6: Obtaining Markers\n")
        subprocess.run(f"python3 {DIR}/bin/03_ObtainMarkers.py -m {args.kmer}",
                       check=True, shell=True)
        # 04_Genotyping.py
        print("\nStep 4/6: Creating Genotype Table\n")
        subprocess.run(f"python3 {DIR}/bin/04_Genotyping.py -m {args.kmer} -x {args.LowCov}",
                       check=True, shell=True)
        # 05_ObtainSegStats.py
        print("\nStep 5/6: Obtaining Segment Statistics\n")
        subprocess.run(f"python3 {DIR}/bin/05_ObtainSegStats.py -m {args.kmer} -L {args.LOD} -f {args.fXX}",
                       check=True, shell=True)

        if (args.Max is not None):
            print("\nExtra Step 5/6: Reducing Number of Markers\n")
            marker_reduction(args.kmer, args.Max)

        # 06_ExportToLepMap3.py
        print("\nStep 6/6: Making Information Usable for LepMap3\n")
        subprocess.run(f"python3 {DIR}/bin/06_ExportToLepMap3.py -m {args.kmer}",
                       check=True, shell=True)
        # 07_LepMap3.py
        print("\nRunning LepMap3\n")
        subprocess.run(f"python3 {DIR}/bin/07_LepMap3.py -m {args.kmer} -t {args.threads} -L {args.LOD} -n {args.nLG}",
                       check=True, shell=True)

        print("AFLAP complete!")

    except subprocess.CalledProcessError as e:
        exit(f"An error occurred: {e}")
