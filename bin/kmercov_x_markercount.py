import pandas as pd
import matplotlib.pyplot as plt

def plot_cov_and_mcount(mc:pd.DataFrame, outfile_name:str)->None:
    x = mc["K-mer Coverage"]
    y = mc["Marker Count"]

    plt.scatter(x, y, c='k')
    plt.xlabel("K-mer Coverage")
    plt.ylabel("Marker Count")

    plt.savefig(outfile_name)
    plt.clf()