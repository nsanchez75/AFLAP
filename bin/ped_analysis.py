import shutil
import os

def sort_ped(f_type:str)->None:
    with open(f"AFLAP_tmp/Pedigree_{f_type}.txt", 'r+') as f:
        lines = sorted(f.readlines())
        f.seek(0)
        f.writelines(lines)

def pedigree_analysis(pedigree: str)->None:
    # try:
        # initialize variables
        parents    = set()
        f1_progs   = set()
        f2_progs   = set()
        f1_crosses = dict()
        f2_crosses = dict()

        # copy pedigree file into AFLAP_Results
        shutil.copy2(pedigree, "AFLAP_Results/Pedigree.txt")        

        # categorize pedigree into F0, F1, F2
        with open(pedigree, 'r') as fin,                        \
             open("AFLAP_tmp/PedigreeInfo.txt", 'w') as fped,   \
             open("AFLAP_tmp/Pedigree_F0.txt", 'w') as f0,      \
             open("AFLAP_tmp/Pedigree_F1.txt", 'w') as f1,      \
             open("AFLAP_tmp/Pedigree_F2.txt", 'w') as f2:

            # write pedigree file into Pedigree.txt
            fped.write(f"Source: {pedigree}")

            # categorize individuals by generation type
            for line in fin:
                cols = line.strip().split()
                # check line size
                if len(cols) != 5:
                    raise ValueError("Line in pedigree file over 5 columns.")

                if cols[1] == '0':
                    if cols[0] not in parents:
                        parents.add(cols[0])

                    f0.write(line)
                elif cols[1] == '1':
                    if cols[0] not in f1_progs:
                        f1_progs.add(cols[0])
                        # create new cross in crosses dictionary
                        cross = f"{cols[3]} {cols[4]}"
                        if cross not in f1_crosses:
                            if f"{cols[4]} {cols[3]}" in f1_crosses:
                                raise ValueError("A parent in the pedigree is being treated as both male and female.")
                            f1_crosses[cross] = 1
                        else: f1_crosses[cross] += 1

                    # check if progeny has duplicated parents
                    if "NA" not in {cols[3], cols[4]} and cols[3] == cols[4]:
                        raise ValueError(f"Identical crossed parents identified for {cols[0]}.")

                    f1.write(line)
                elif cols[1] == '2':
                    if cols[0] not in f2_progs:
                        f2_progs.add(cols[0])

                    f2.write(line)
                else:
                    raise ValueError("Pedigree file contains individual that is not F0, F1, or F2.")


        # check Pedigree_F1.txt and Pedigree_F2.txt structure
        if os.path.exists("AFLAP_tmp/Pedigree_F1.txt"):
            with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as f:
                for line in f:
                    line = line.strip().split()
                    if line[3] not in parents or line[4] not in parents:
                        raise ValueError(f"F1 progeny {line[0]} descends from a parent not found in pedigree file.")
        if os.path.exists("AFLAP_tmp/Pedigree_F2.txt"):
            with open("AFLAP_tmp/Pedigree_F2.txt", 'r') as f:
                for line in f:
                    line = line.strip().split()
                    if line[3] not in f1_progs or line[4] not in f1_progs:
                        raise ValueError(f"F2 progeny {line[0]} descends from an F1 progeny not found in pedigree file.")

        # sort pedigree files
        sort_ped("F0")
        sort_ped("F1")
        sort_ped("F2")

        # determine parent count
        if len(parents) < 2:
            raise ValueError("Pedigree file does not have at least 2 parents.")
        elif len(parents) == 2: print("2 parents detected. This will be easy! Identifying cross(es)...")
        else: print(f"{len(parents)} parents detected. This will not be so easy! Identifying crosses...")

        # identify and print cross
        with open("AFLAP_tmp/Crosses.txt", 'w') as f:
            ## F1 progeny
            print("F1 crosses that have been identified:")
            for cross in f1_crosses:
                c_vals = cross.strip().split()
                print(f"\t{c_vals[0]}x{c_vals[1]}")
                f.write(f"{f1_crosses[cross]} 1 {c_vals[0]} {c_vals[1]}")
            print()
            ## F2 progeny
            # TODO: implement stuff below when working with F2
            print("F2 crosses that have been identified:")
            for crosses in f2_crosses:
                print(f"\t{crosses[0]}x{crosses[1]}")
                f.write(f"{f2_crosses[crosses]} 2 {crosses[0]} {crosses[1]}")

        # perform LA analysis
        checked_parents = set()
        with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as fin, \
             open("AFLAP_tmp/LA.txt", 'w') as fla,       \
             open("AFLAP_tmp/noLA.txt", 'w') as fnola:

            for line in fin:
                line = line.strip().split()
                if line[0] not in checked_parents:
                    if "NA" in (line[3], line[4]):
                        fnola.write(f"{line[0]}\n")
                        checked_parents.add(line[0])
                    elif [isinstance(x, int) for x in [line[3], line[4]]]:
                        if line[3] > line[4]:
                            raise("Cannot have a lower bound higher than an upper bound.")

                        # write {parent} {lower bound} {upper bound} to LA.txt
                        fla.write(f"{line[0]} {line[3]} {line[4]}")
                        checked_parents.add(line[0])
                    else:
                        raise ValueError("Invalid bound entry in Pedigree_F0.txt.")

        print("Finished pedigree file analysis.")

    # except Exception as e:
    #     print(f"Error when running ped_analysis: {e}")
    #     exit(1)