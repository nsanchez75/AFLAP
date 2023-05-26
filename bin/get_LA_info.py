import os

def get_LA_info()->list[tuple[str, int, int, str]]:
    la_filename = "AFLAP_tmp/LA.txt"
    cross_filename = "AFLAP_tmp/Crosses.txt"

    # check if files exist
    for file in (la_filename, cross_filename):
        if (not os.path.exists(file)):
            raise FileNotFoundError(f"Error: {file} not found.")

    with open(la_filename, 'r') as fla:
        ret_list = list()
        for p in fla:
            p = p.strip().split()
            G = p[0]
            LO = int(p[1])
            UP = int(p[2])

            p0 = list()
            with open(cross_filename, 'r') as fc:
                for pair in fc:
                    pair = pair.strip().split()

                    if pair[2] == G:
                        p0.append(pair[3])
                    elif pair[3] == G:
                        p0.append(pair[2])
                    else: continue

            p0 = '_'.join(p0)

            # create tuple
            ret_list.append((G, LO, UP, p0))

    return ret_list