from numpy import*

def readparams(filename):
    """
    Read tipsy ".param" file and return a dictionary 
    """

    #read each line in file and ignor blank ones
    with open(filename) as FILE:
        lines = list(filter(None,(line.rstrip() for line in FILE))) 
    
    #store name and values, ignore commented lines and spaces here
    names, values = [line.split('=')[0] for line in lines if line[0]!='#'], [line.split('=')[1] for line in lines if line[0]!='#']
    names         = [n.strip() for n in names]
    values        = [v.strip() for v in values]

    #params is str:str for now...
    params        = dict(zip(names, values))

    return params
