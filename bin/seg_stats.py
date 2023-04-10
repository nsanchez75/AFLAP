import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def normalize_data(markers:pd.DataFrame)->pd.DataFrame:
    return markers["Frequency Count"] / markers["Frequency Count"].sum()

def get_datapoints(markers:pd.DataFrame)->tuple[np.ndarray, np.ndarray]:
    return markers["Frequency"].to_numpy(), markers["Frequency Count"].to_numpy()

def get_seg_stats(markers_all:pd.DataFrame, markers_equal:pd.DataFrame, markers_over:pd.DataFrame, ak:int, oufile_name:str)->None:
    # normalize data
    markers_all["Frequency Count"]      = normalize_data(markers_all)
    markers_equal["Frequency Count"]    = normalize_data(markers_equal)
    markers_over["Frequency Count"]     = normalize_data(markers_over)

    x1, y1 = get_datapoints(markers_equal)
    x2, y2 = get_datapoints(markers_over)
    x3, y3 = get_datapoints(markers_all)

    # init plot
    fig, ax = plt.subplots()

    # plot data
    ax.plot(x1, y1, 'ro', label=f'={ak}')
    ax.plot(x2, y2,'bo', label=f'>{ak}')
    ax.plot(x3, y3, 'k-')

    # implement plot info
    plt.xlabel("Marker Presence")
    plt.ylabel("Marker Density")
    ax.legend('best')

    # create png
    plt.savefig(oufile_name)