import matplotlib.pyplot as plt
import argparse
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

parser = argparse.ArgumentParser(description='Routine to read and plot disperse persistence pairs. Args: < NDnet file ascii format (str) >, < optional -- second hdf5 file (str)>, < suffix of outputfile (str) >, < what to plot (str) >, < plane to plot (str: XY, YZ, XZ)>, < width of slice (int, in cells)>, <optional -- max T (float, K)>, <optional -- min T (float, K)>')
parser.add_argument('fname',nargs=1,type=str, help='Absolute path to NDnet file + filename (ascii)')
parser.add_argument('--cut',nargs=1,type=bool, help='Plot lines at given cut')
parser.add_argument('--suffix',nargs=1,type=str, help='Suffix of outputfile (str, optional)')
parser.add_argument('toplot',nargs=1,type=str, help='Quantity to plot')
parser.add_argument('plane',nargs=1,type=str, help='Plane to plot')
parser.add_argument('width',nargs=1,type=int, help='Width of slice (grid cells)')
parser.add_argument('--tempmin',nargs=1,type=float, help='Minimum temperature (float, optional)')
parser.add_argument('--tempmax',nargs=1,type=float, help='Maximum temperature (float, optional)')
args   = parser.parse_args()
path   = args.hdf5file[0]

def read_ppairs_ascii(fname, cut = False, sigma= True, prob_x = 0.999, prob_y = 0.999, plot_prob_x = True, plot_prob_y = True, sigma_x = 5., sigma_y=5. ):
    with open(fname, 'r') as file:
        lines = file.read().splitlines()
    
    sections = ['type', 'persistence', 'field_value']

    idxs = []
    for i in range(len(sections)):
        idxs.append(lines.index(sections[i]))
    
    n_points = int(lines[4])
    print('# Persistence Pairs: ',n_points)

    type = lines[idxs[0]+1 : idxs[0] + n_points +1]
    color = [float(p) for p in type]

    persistence = lines[idxs[1]+1 : idxs[1] + n_points +1]
    y = [float(p) for p in persistence]

    field = lines[idxs[2]+1 : idxs[2] + n_points+1]
    x = [float(f) for f in field]

    type = unique(color)

    y_3 = [y[i] for i in range(len(y)) if color[i]==3.]
    x_1 = [x[i] for i in range(len(x)) if color[i]==1.]

    if cut:
        quantile_x_1 = quantile(x_1, prob_x) #0.999 \approx 3.5 sigma
        quantile_y_3 = quantile(y_3, prob_y) #0.9995 \approx 3.5 sigma
    elif sigma:
        quantile_x_1 = std(x)
        quantile_x_1 = sigma_x * quantile_x_1

        quantile_y_3 = std(y)
        quantile_y_3 = sigma_y * quantile_y_3


    #--Plot--#
    fig, ax = plt.subplots()

    scatter = ax.scatter(x, y, c=color, s=0.3, alpha = 0.1)

    legend = ax.legend(*scatter.legend_elements(),
                    loc="best", title="Critical Points")
    for i in legend.legendHandles:
        i.set_alpha(1)
    ax.add_artist(legend)

    if plot_prob_y:
        plt.axhline(quantile_y_3, min(x), max(x), color = 'black')
    if plot_prob_x:
        plt.axvline(quantile_x_1, min(y), max(y), color='black')
    plt.xlabel('Field')
    plt.ylabel('Persistence')
    plt.loglog()

    fig.savefig(c.PLOTPATH + '/ppairs.pdf')

    return quantile_x_1, quantile_y_3
