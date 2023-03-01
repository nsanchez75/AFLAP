import os


def check_01_jellyfish():
    if not os.path.exists("AFLAP_tmp/01/LA.txt"):
        print("Could not find AFLAP_tmp/01/LA.txt. Please rerun 01_JELLYFISH.py. Terminating.")
        exit(1)

def check_02_escm():
    if not os.path.exists("AFLAP_tmp/02/Boundaries.txt"):
        print("Could not find AFLAP_tmp/02/Boundaries.txt. Please rerun 02_ExtractSingleCopyMers.py. Terminating.")
        exit(1)
    