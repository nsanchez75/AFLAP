import os

def sort_ped(f_type:str)->None:
    with open(f"AFLAP_tmp/Pedigree_{f_type}.txt", 'r+') as f:
        lines = sorted(f.readlines())
        f.seek(0)
        f.writelines(lines)

def sortwrite_Ftext(in_set:set, f_type:str)->None:
    with open(f"AFLAP_tmp/01/{f_type}.txt", 'r+') as f:
        for s in in_set: f.write(f"{s}\n")
        f.seek(0)
        lines = sorted(f.readlines())
        f.seek(0)
        f.writelines(lines)

def pedigree_analysis(pedigree: str)->None:
    # 0. initialize variables
    f0_count   = 0
    f1_count   = 0
    f2_count   = 0
    parents    = set()
    f1_progs   = set()
    f2_progs   = set()
    f1_crosses = dict()
    f2_crosses = dict()


    # 1. categorize pedigree into F0, F1, F2
    with open(pedigree, 'r') as fin, open("AFLAP_tmp/Pedigree.txt", 'w') as fped, open("AFLAP_tmp/Pedigree_F0.txt", 'w') as f0, open("AFLAP_tmp/Pedigree_F1.txt", 'w') as f1, open("AFLAP_tmp/Pedigree_F2.txt", 'w') as f2:
        # write pedigree file into Pedigree.txt
        fped.write(f"Original File Source: {pedigree}")

        for line in fin:
            cols = line.strip().split()
            # check line size
            if len(cols) != 5:
                raise ValueError("Error: Line in pedigree file over 5 columns.")

            if cols[1] == '0':
                if cols[0] not in parents:
                    f0_count += 1
                    parents.add(cols[0])

                f0.write(line)
            elif cols[1] == '1':
                if cols[0] not in f1_progs:
                    f1_count += 1
                    f1_progs.add(cols[0])
                    # add or create new cross in crosses dictionary
                    cross = f"{cols[3]} {cols[4]}"
                    if cross not in f1_crosses:
                        if f"{cols[4]} {cols[3]}" in f1_crosses:
                            raise ValueError("Error: A parent in the pedigree is being treated as both male and female.")
                        f1_crosses[cross] = 1
                    else: f1_crosses[cross] += 1

                # check if progeny has duplicated parents
                if "NA" not in {cols[3], cols[4]} and cols[3] == cols[4]:
                    raise ValueError(f"Error: Identical crossed parents identified for {cols[0]}.")

                f1.write(line)
            elif cols[1] == '2':
                # TODO: work on later once we begin working on F2
                raise ValueError("Error: AFLAP does not currently work on F2 progeny.")
                f2.write(line) # for future use
            else:
                raise ValueError("Error: Pedigree file contains individual that is not F0, F1, or F2.")


    # 2. check F1.txt and F2.txt structure
    if os.path.exists("AFLAP_tmp/Pedigree_F1.txt"):
        with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as f:
            for line in f:
                line = line.strip().split()
                if line[3] not in parents or line[4] not in parents:
                    raise ValueError(f"Error: F1 progeny {line[0]} descends from parent not found in pedigree file.")
    if os.path.exists("AFLAP_tmp/Pedigree_F2.txt"):
        pass    # TODO: implement when working w/ F2


    # 3. sort pedigree files
    sort_ped("F0")
    sort_ped("F1")
    sort_ped("F2")


    # 4. create F0.txt, F1.txt, F2.txt  FIXME: could be redundant with Pedigree_{f_type}.txt having same info
    sortwrite_Ftext(parents, "F0")
    sortwrite_Ftext(f1_progs, "F1")
    # sortwrite_Ftext(f2_progs, "F2")   # TODO: implement F2.txt when working w/ F2


    # 5. determine parent count
    if f0_count < 2:
        raise ValueError("Error: Pedigree file does not identify at least 2 parents.")
    elif f0_count == 2: print("2 parents detected. This will be easy! Identifying cross(es)...")
    else: print(f"{f0_count} parents detected. This will not be so easy! Identifying crosses...")


    # 6. identify and print cross
    with open("AFLAP_tmp/01/Crosses.txt", 'w') as f:
        # F1 progeny
        print("F1 crosses that have been identified:")
        for cross in f1_crosses:
            c_vals = cross.strip().split()
            print(f"\t{c_vals[0]}x{c_vals[1]}")
            f.write(f"{f1_crosses[cross]} 1 {c_vals[0]} {c_vals[1]}")
        print()
        # F2 progeny
        # # TODO: implement stuff below when working with F2
        # print("F1 crosses that have been identified:")
        # for crosses in f2_crosses:
        #     print(f"\t{crosses[0]}x{crosses[1]}")
        #     f.write(f"{f2_crosses[crosses]} 2 {crosses[0]} {crosses[1]}")