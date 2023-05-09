from matplotlib import pyplot as plt
import os

def find_max_y(yvals:list[int], LO:int, HI:int)->int:
    max_y = 0
    if len(yvals) < LO: max_y = max(yvals)
    else:
        max_y = 0
        for i in range(LO, HI):
            if max_y < yvals[i]: max_y = yvals[i]

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
    max_y = find_max_y(histo_y, LO, HI)

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