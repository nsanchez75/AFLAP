import subprocess

def separate_ped(pedigree):
    with open(pedigree, 'r') as fin, open("AFLAP_tmp/Pedigree_F0.txt", 'w') as f0, open("AFLAP_tmp/Pedigree_F1.txt", 'w') as f1, open("AFLAP_tmp/Pedigree_F1.txt", 'w') as f2:
        labels = dict()
        for line in fin:
            col = line.strip().split()

            # quick label duplicate check
            if col[0] in labels and labels[col[0]] != col[1]:
                print("Pedigree file specifies the same individual label over multiple generations. " +
                      "This is incompatable with AFLAP. Terminating.")
                exit(1)
            else:
                labels[col[0]] = col[1]

            # separate individuals by generation
            if col[1] == '0':   f0.write(line)
            elif col[1] == '1': f1.write(line)
            elif col[1] == '2': f2.write(line)
            else:
                print("Pedigree file specifies individuals which are not F0, F1 or F2. " +
                      "AFLAP has not been validated on this type of data. Terminating.")
                exit(1)


def check_format(filename)->None:
    with open(filename, 'r') as f:
        for line in f:
            col = line.strip().split()
            if len(col) != 5:
                print("Pedigree file not formatted correctly! Five columns should be present on " +
                      "every progeny line with the fourth and fifth specifying the parents. " +
                      "E.g:\n\t[ProgLabel]\t[1/2]\t[READS]\t[Parent1]\t[Parent2]\n" +
                      "Terminating.")
                exit(1)


def confirm_parents()->None:
    with open("AFLAP_tmp/Pedigree_F0.txt", 'r') as fin, open("AFLAP_tmp/01/F0.txt", 'w+') as fout:
        parents = set()
        for line in fin:
            col = line.strip().split()
            # find parents (0th generation) in file
            if col[0] not in parents:
                fout.write(f"{col[0]}\n")
                parents.add(col[0])

        fout.seek(0)
        p_count = len(fout.readlines())
        if p_count < 2:
            print(f"{p_count} parents detected. Two or more F0 should be specified in the pedigree file to run AFLAP. Terminating.")
            exit(1)
        elif p_count == 2:
            print(f"{p_count} parents detected. Simple! Beginning k-mer counting...")
        else:
            print(f"{p_count} parents detected. Not so simple! Determining crosses...")

    with open("AFLAP_tmp/Pedigree_F1.txt", 'r') as fin, open("AFLAP_tmp/01/Crosses.txt", 'w') as fout:
        crosses = dict()
        labels = set()
        for line in fin:
            col = line.strip().split()
            cross = f"{col[3]} {col[4]}"

            # FIXME: add check for if there are two different parents for same progeny (shouldn't be possible)
            if col[0] not in labels:
                if cross not in crosses:
                    crosses[cross] = 1
                else:
                    crosses[cross] += 1
                    labels.add(col[0])
                labels.add(col[0])

        for i in crosses:
            fout.write(f"{crosses[i]} 1 {i}\n")
    try:
        with open("AFLAP_tmp/01/Crosses.txt", 'r') as f:
            crosses = f.readlines()
            cross_count = len(crosses)
            if cross_count > 4:
                print(f"{cross_count} crosses identified in pedigree file. This may take a long time to process. Crosses include:")
            else:
                print(f"{cross_count} cross(es) identified:")
            
            for cross in crosses:
                cross = cross.strip().split()
                print(f"\t\t{cross[2]} {cross[3]}")
    except FileNotFoundError:
        print("Could not find AFLAP_tmp/01/Crosses.txt. Terminating.")
        exit(1)


def check_pedigree_structure()->None:
    with open("AFLAP_tmp/01/Crosses.txt", 'r') as f:
        f1 = 0
        f2 = 0

        for line in f:
            col = line.strip().split()
            if int(col[1]) == 1:
                f1 += 1
            elif int(col[1]) == 2:
                f2 += 1

            if f1 and f2:
                print("Both F1 and F2 populations detected. AFLAP has not been validated on this type of data so will terminate. " +
                      "Please try to rerun specifying either type of population. Terminating.")
                exit(1)

        if not f1 and not f2:
            print("Neither F1 or F2 populations detected. AFLAP requires either F1 or F2 populations in order to run. Terminating.")
            exit(1)
        else:
            if f1 and not f2:
                print(f"\t{f1} F1 cross(es) identified.")
                temp_val = 1

            else: # f2 no f1
                print(f"\t{f2} F2 cross(es) identified")
                temp_val = 2

            with open("AFLAP_tmp/01/Crosses.txt", 'r') as fin, open("AFLAP_tmp/01/ParentsToCompare.txt", 'w+') as fout:
                for line in fin:
                    columns = line.strip().split()
                    if int(columns[1]) == temp_val:
                        fout.write(f"{columns[2]}x{columns[3]}\n")