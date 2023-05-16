import argparse
import os
import pandas as pd
import multiprocessing as mp
import subprocess

import get_LA_info as gli

def run_subprocess(cmd, outfile, errfile):
    outf = open(outfile, 'w')
    errf = open(errfile, 'w')
    subprocess.Popen(cmd, shell=True, stdout=outf, stderr=errf)
    outf.close()
    errf.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='LepMap3', description='A script to run LepMap3 and produce a genetic map which can be aligned to a genome assembly.')
    parser.add_argument('-m', '--kmer', type=int, default=31, help='K-mer size (optional). Default [31].')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for JELLYFISH counting (optional). Default [4].')
    parser.add_argument('-L', '--LOD', type=int, default=2, help='LOD score - Will run LepMap3 with minimum LOD. Default [2].')
    args = parser.parse_args()

    # check if necessary files exist
    if not os.path.exists("AFLAP_tmp/01/LA.txt"):
        raise FileNotFoundError("Error: AFLAP_tmp/01/LA.txt not found.")

    # TODO: check if $CONDA_PREFIX exists

    # make directory
    os.makedirs(f"AFLAP_Results/LOD{args.LOD}", exist_ok=True)

    # run LepMap3 on parents
    list_of_Gs = gli.get_LA_info()
    sc2_processes = list()
    for G_info in list_of_Gs:
        G, LO, UP, P0 = G_info

        # check if tsv file for G exists
        if not os.path.exists(f"AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv"):
            raise FileNotFoundError(f"Error: LepMap-ready genotype table fo {G} not found.")

        # begin LepMap3 analysis
        print(f"Initiating LepMap3 analysis on {G}...")
        ## skip if G had already been analyzed
        if os.path.exists(f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.txt"):
            print(f"\tPrevious results for {G} detected. Skipping.")
        else:
            # run LepMap3 - SeparateChromosomes2
            sc2_stderr = f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.stderr"
            sc2_stdout = f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.txt"
            p = mp.Process(target=run_subprocess, args=(f"java -cp $CONDA_PREFIX/bin/lepmap3/ SeparateChromosomes2 lodLimit={args.LOD} numThreads={args.threads} data=AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv", sc2_stdout, sc2_stderr)).start()
            sc2_processes.append(p)
            # subprocess.run(args=f"java -cp $CONDA_PREFIX/bin/lepmap3/ SeparateChromosomes2 lodLimit={args.LOD} numThreads={args.threads} data=AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv",
            #                stdout=sc2_stdout, stderr=sc2_stderr, shell=True)

        # wait for SeparateChromosomes2 processes to finish
        for p in sc2_processes:
            mp.Process(p).join()
        print(f"Finished analyzing {G}.")


        # gather analysis statistics
        m_count = 0
        fre_dict = dict()
        with open(f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.txt", 'r') as flod:
            for res_line in flod:
                # skip file header
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

        if (args.threads > len(lg_set)):
            print("Running linkage group ordering in parallel as number of threads exceeds number of linkage groups...")
        else:
            print("Running linkage group ordering in series as number of threads does not exceed number of linkage groups...")
        
        om2_processes = list()
        for lg in lg_set:
            if os.path.exists(f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.LG{lg}.txt"):
                print(f"\tAnalysis of linkage group {lg} detected. Skipping.")
            else:
                om2_stdout = f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.LG{lg}.txt"
                om2_stderr = f"AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.LG{lg}.stderr"
                p = mp.Process(target=run_subprocess, args=(f"java -cp $CONDA_PREFIX/bin/lepmap3/ OrderMarkers2 useMorgan=1 numMergeIterations=20 chromosome={lg} map=AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.txt data=AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv", om2_stdout, om2_stderr)).start()
                # subprocess.run(args=f"java -cp $CONDA_PREFIX/bin/lepmap3/ OrderMarkers2 useMorgan=1 numMergeIterations=20 chromosome={lg} map=AFLAP_Results/LOD{args.LOD}/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.LOD{args.LOD}.txt data=AFLAP_Results/{G}_m{args.kmer}_L{LO}_U{UP}_{P0}.ForLepMap3.tsv",
                #                stdout=om2_stdout, shell=True)

                print(f"\tAnalysis of linkage group {lg} complete.")

        print("Linkage group ordering complete")

        # with open("AFLAP_tmp/01/Crosses.txt", 'r') as fcrosses:
        #     for cross in fcrosses:
        #         cross = cross.strip().split()
        #         if (cross[2] == G):
        #             sex_check = 0
        #         elif (cross[3] == G):
        #             sex_check = 1
        #         else:
        #             continue        
        #         break
        
        # if (not sex_check):


        print("continue coding 07?")