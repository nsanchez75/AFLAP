import argparse
import glob
import pandas as pd
import os
import multiprocessing as mp
import subprocess

from get_LA_info import get_LA_info

#################################################
#       A Python script to run LepMap3 and produce a genetic map which can be aligned to a genome assembly.
#################################################

def run_sc2(lodfile:str, loderr:str, forlepmap_file:str, LOD:int, parent_header:str, num_threads:int)->None:
    if os.path.exists(lodfile):
        print(f"\tPrevious {LOD}-score results for {parent_header} detected. Skipping.")
        return

    with open(lodfile, 'w') as sc2_stdout, open(loderr, 'w') as sc2_stderr:
        subprocess.run(args=f"java -cp $CONDA_PREFIX/bin/lepmap3/ SeparateChromosomes2 lodLimit={LOD} numThreads={num_threads} data={forlepmap_file}",
                       stdout=sc2_stdout, stderr=sc2_stderr, shell=True)

def run_om2(txt_header:str, LOD:int, lg:int, outfile:str, errfile:str)->None:
    with open(outfile, 'w') as outf, open(errfile, 'w') as errf:
        subprocess.run(args=f"java -cp $CONDA_PREFIX/bin/lepmap3/ OrderMarkers2 useMorgan=1 numMergeIterations=20 chromosome={lg} map=AFLAP_Results/LOD{LOD}/{txt_header}.LOD{LOD}.txt data=AFLAP_Results/{txt_header}.ForLepMap3.tsv",
                       stdout=outf, stderr=errf, shell=True)
    print(f"\tAnalysis of linkage group {lg} complete.")

def run_lepmap(parent_header:str, txt_header:str, f_type:str, num_threads:int, LOD:int, SEX:str=None, num_LGs:int=None)->None:
    # check if tsv file for G exists
    forlepmap_file = f"AFLAP_Results/{txt_header}.ForLepMap3.tsv"
    if not os.path.exists(forlepmap_file):
        exit(f"An error occurred: LepMap-ready genotype table for {parent_header} not found. Rerun 06_ExportToLepMap3.py.")

    # begin LepMap3 analysis
    print(f"Initiating LepMap3 analysis on {f_type} progeny of {parent_header}...")
    print(f"\tWarning: If running an analysis on F2 populations, note that the program will run continuously until at least {num_LGs} linkage groups are found.")
    while True:
        os.makedirs(f"AFLAP_Results/LOD{LOD}/{f_type}", exist_ok=True)
        lodfile = f"AFLAP_Results/LOD{LOD}/{f_type}/{txt_header}.LOD{LOD}.txt"
        loderr = f"AFLAP_Results/LOD{LOD}/{f_type}/{txt_header}.LOD{LOD}.stderr"
        # run LepMap3 - SeparateChromosomes2
        run_sc2(lodfile, loderr, forlepmap_file, LOD, parent_header, num_threads)
        with open(lodfile, 'r') as flod: sc2_results = set(int(line.strip()) for line in flod.readlines() if not (line.strip().startswith('#')))
        # stop analysis if requirements satisfied
        if f_type == "F1" or len(sc2_results) >= num_LGs:
            break
        LOD += 1
    print(f"{len(sc2_results)} linkage groups found.")

    # gather analysis statistics
    m_count = 0
    fre_dict = dict()
    with open(lodfile, 'r') as flod:
        for res_line in flod:
            if res_line.startswith('#'): continue
            res_line = int(res_line.strip())
            if res_line not in fre_dict: fre_dict[res_line] = 1
            else: fre_dict[res_line] += 1

            m_count += 1
    ## count number of linkage groups
    lg_set = set()
    for res in fre_dict:
        if res and (fre_dict[res]/m_count >= 0.01):
            lg_set.add(res)
    print(f"{len(lg_set)} linkage groups detected containing a minimum of 1% of the markers.")

    # analyze threads
    print("Running linkage group ordering in", sep=' ')
    if (num_threads > len(lg_set)): print("parallel as number of threads exceeds number of linkage groups...")
    else: print("series as number of threads does not exceed number of linkage groups...")

    # process linkage groups
    om2_processes = list()
    for lg in lg_set:
        if os.path.exists(f"AFLAP_Results/LOD{LOD}/{f_type}/{txt_header}.LOD{LOD}.LG{lg}.txt"):
            print(f"\tAnalysis of linkage group {lg} detected. Skipping.")
        else:
            om2_stdout = f"AFLAP_Results/LOD{LOD}/{f_type}/{txt_header}.LOD{LOD}.LG{lg}.txt"
            om2_stderr = f"AFLAP_Results/LOD{LOD}/{f_type}/{txt_header}.LOD{LOD}.LG{lg}.stderr"
            p = mp.Process(target=run_om2, args=(txt_header, LOD, lg, om2_stdout, om2_stderr))
            p.start()
            om2_processes.append(p)
    for p in om2_processes:
        p.join()
    print("Linkage group ordering complete.")

    # get row indices of marker sequences (used filtered genotype table as reference)
    filtered_tsv = f"AFLAP_tmp/05/{txt_header}.Genotypes.MarkerID.Filtered.tsv"
    if not os.path.exists(filtered_tsv):
        exit(f"An error occurred: {filtered_tsv} not found. Rerun 05_ExportToLepMap3.py")
    markerid_df = pd.read_csv(filtered_tsv, sep='\t', usecols=["MarkerSequence", "MarkerID"])
    markerid_df["RowIndex"] = markerid_df.index + 1

    # identify linkage groups to marker sequences
    lg_df = pd.DataFrame()
    for glob_path in glob.glob(f"AFLAP_Results/LOD{LOD}/{f_type}/{txt_header}.LOD{LOD}.LG*.txt"):
        if SEX == "male":     COLS_USED, COL_NAMES = [0, 1], ["RowIndex", "Position"]
        elif SEX == "female": COLS_USED, COL_NAMES = [0, 2], ["RowIndex", "Position"]
        ## note: above 2 conditions may not pass since SEX F1-specific
        else:                 COLS_USED, COL_NAMES = [0, 1, 2], ["RowIndex", "Male Position", "Female Position"]
        lg_info = pd.read_csv(glob_path, sep='\t', names=COL_NAMES, skiprows=3, usecols=COLS_USED)
        lg_info['LG'] = glob_path.removeprefix(f"AFLAP_Results/LOD{LOD}/{f_type}/{txt_header}.LOD{LOD}.LG").removesuffix(".txt")

        lg_df = lg_info if (lg_df.empty) else pd.concat([lg_df, lg_info])

    joined_df = pd.merge(markerid_df, lg_df, on="RowIndex", how='inner')
    joined_df.to_csv(f"AFLAP_Results/{txt_header}.LOD{LOD}.txt", sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='LepMap3', description='A script to run LepMap3 and produce a genetic map which can be aligned to a genome assembly.')
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size. Default [31].')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for JELLYFISH counting. Default [4].')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    parser.add_argument('-n', '--nLG', type=int, default=10, help='Minimum number of linkage groups needed to continue F2 LepMap3 analysis. Default [10].')
    args = parser.parse_args()

    if not os.environ.get("CONDA_PREFIX"):
        exit("An error occurred: $CONDA_PREFIX does not exist. Activate conda and ensure the AFLAP.yml environment is satisfied.")

    for f_type in ["F1", "F2"]:
        if not len(os.listdir(f"AFLAP_tmp/01/{f_type}Count")):
            print(f"No {f_type} progeny found. Skipping.")
            continue

        if f_type == "F1":
            for G, LO, UP, P0, SEX in get_LA_info():
                parent_header = G
                txt_header = f"{G}_F1_m{args.kmer}_L{LO}_U{UP}_{P0}"
                run_lepmap(parent_header, txt_header, f_type, args.threads, args.LOD, SEX=SEX)
        else:
            male = list()
            female = list()
            for G, _, _, _, SEX in get_LA_info():
                if SEX == "male": male.append(G)
                else:             female.append(G)
            male = '_'.join(male)
            female = '_'.join(female)

            parent_header = f"{male}x{female}"
            txt_header = f"{parent_header}_F2_m{args.kmer}"
            run_lepmap(parent_header, txt_header, f_type, args.threads, args.LOD, num_LGs=args.nLG)
