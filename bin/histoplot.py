from matplotlib import pyplot as plt
import os

def histoplot(infilepath:str, LO:int, HI:int, outfilepath:str)->None:
    # check if file exists
    if not os.path.exists(infilepath):
        raise FileNotFoundError(f"Error: {infilepath} not found.")

    histo_x = list()
    histo_y = list()

    # get data
    with open(infilepath) as fin:
        for line in fin:
            line = line.strip().split()
            histo_x.append(int(line[0]))
            histo_y.append(int(line[1]))

    # find maximum y bound
    max_y = 0
    print(f"Length of histo_y: {len(histo_y)}")
    print(f"LO: {LO}\nHI: {HI}")
    if len(histo_y) < LO: max_y = max(histo_y)
    else:
        max_y = 0
        for i in range(LO, HI):
            if max_y < histo_y[i]: max_y = histo_y[i]

    # plot histogram
    plt.plot(histo_x, histo_y, 'k')
    plt.vlines(x=[LO, HI], ymin=0, ymax=2*max_y, linestyles='dashed', colors='k')
    plt.xlim(0, 300)
    plt.ylim(0, 2*max_y)
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")

    # create png
    plt.savefig(outfilepath)
    plt.clf()