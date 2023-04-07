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
    markers_all     = normalize_data(markers_all)
    markers_equal   = normalize_data(markers_equal)
    markers_over    = normalize_data(markers_over)

    x1, y1 = get_datapoints(markers_equal)
    x2, y2 = get_datapoints(markers_over)

    