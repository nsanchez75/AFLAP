import pandas as pd
import matplotlib.pyplot as plt
import sys

def normalize_data(markers:pd.DataFrame)->pd.DataFrame:
    return markers["Frequency Count"] / markers["Frequency Count"].sum()

def get_seg_stats(markers_all:pd.DataFrame, markers_equal:pd.DataFrame, markers_over:pd.DataFrame, oufile_name:str)->None:
    # normalize data
    markers_all     = normalize_data(markers_all)
    markers_equal   = normalize_data(markers_equal)
    markers_over    = normalize_data(markers_over)

    x1 = markers_equal["Frequency"].to_numpy()