import numpy as np
import matplotlib.pyplot as plt
#import constants as c

def read_ppairs_ascii(fname, cut = False, sigma= True, prob_x = 0.999, prob_y = 0.999, plot_prob_x = True, plot_prob_y = True, sigma_x = 5., sigma_y=5. ):
    with open(fname, 'r') as file:
        lines = file.read().splitlines()
    
    sections = ['type', 'persistence', 'field_value']

    idxs = []
    for i in range(len(sections)):
        idxs.append(lines.index(sections[i]))
   
    print(idxs)

    n_points = int(lines[4])
    print('# Persistence Pairs: ',n_points)

    type = lines[idxs[0]+1 : idxs[0] + n_points +1]
    subidxs = [t for t in range(len(type)) if type[t] == '2']
    n_points_2 = len(subidxs)

    
    type = lines[idxs[0]+1 : idxs[0] + n_points_2 +1]
    color = [float(p) for p in type]

    persistence = lines[idxs[1]+1 : idxs[1] + n_points_2 +1]
    y = [float(p) for p in persistence]

    field = lines[idxs[2]+1 : idxs[2] + n_points_2+1]
    x = [float(f) for f in field]

    type = np.unique(color)

    y_3 = [y[i] for i in range(len(y)) if color[i]==3.]
    x_1 = [x[i] for i in range(len(x)) if color[i]==1.]

    if cut:
        quantile_x_1 = np.quantile(x_1, prob_x) #0.999 \approx 3.5 sigma
        quantile_y_3 = np.quantile(y_3, prob_y) #0.9995 \approx 3.5 sigma
    elif sigma:
        quantile_x_1 = np.std(x)
        quantile_x_1 = sigma_x * quantile_x_1

        quantile_y_3 = np.std(y)
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

    fig.savefig('ppairs.pdf')

    return quantile_x_1, quantile_y_3
