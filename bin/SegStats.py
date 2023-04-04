import pandas as pd
import matplotlib
import sys

if __name__ == "__main__":
    M61 = pd.read_table(sys.argv[1], sep=' ', header=False)
    M62 = pd.read_table(sys.argv[2], sep=' ', header=False)
    Mal1 = pd.read_table(sys.argv[3], sep=' ', header=False)

    M61sum = M61.sum(axis='columns')