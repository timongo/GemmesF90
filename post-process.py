"""
This script can be used to post process the *.out file.
"""
# imports ---------------------------------------------------------------------
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys, os

# functions -------------------------------------------------------------------
def load_args():
    """
    Function that loads data
    """
    argl = sys.argv

    # try to find help
    for arg in argl:
        if arg=='-h':
            fusage()

    # find arguments
    argdict = ARGD
    for n, arg in enumerate(argl):
        if arg[0]=='-':
            Arg = arg[1:]
            if Arg=='name':
                try:
                    argdict[Arg] = argl[n+1]
                except IndexError:
                    ferror(key=1, msg=Arg)
            elif Arg=='plot':
                argdict['plot'] = True
            elif Arg=='doption':
                try:
                    tmp = argl[n+1]
                    varfound = False 
                    for i, var in enumerate(VARS):
                        if tmp==var:
                            varfound = True
                            argdict['col'] = i
                            break
                    argdict['doption'] = tmp if varfound else ferror(key=2)
                except IndexError:
                    ferror(key=1, msg=Arg)

    print('\n> Arguments loaded')
    return argdict

def load_data():
    """
    Fonction to load data from gemmes.
    """
    data = np.loadtxt(fname='gemmes.out', skiprows=1)

    print('\n> Gemmes data loaded')
    return data

def fusage():
    """
    Function Usage
    """
    print('\n> Usage: python post-process.py [-name name] [-doption option]')
    print('> [-h] [-plot]')
    print('>\n> Options:')
    print('>     -h :: display this message.')
    print('>     -name : string :: name of figure.')
    print('>     -plot :: to show the figures instead to save them.')
    print('>     -doption : string :: name of the considerd option in:', VARS)
    print()
    ferror()

def change_mpl_opt(opt):
    """
    This function changes matplotlib options.
    """
    size = SIZE
    figsize = FIGSIZE
    if opt!='all':
        size = 12
        figsize = (5, 5)

    SMALL_SIZE = SIZE 
    MEDIUM_SIZE = SIZE 
    BIGGER_SIZE = SIZE 
    plt.rc('font', size = SMALL_SIZE) # controls default text sizes
    plt.rc('axes', titlesize = BIGGER_SIZE) # fontsize of the axes title
    plt.rc('axes', labelsize = MEDIUM_SIZE) # fontsize of the x and y labels
    plt.rc('xtick', labelsize = SMALL_SIZE) # fontsize of the tick labels
    plt.rc('ytick', labelsize = SMALL_SIZE) # fontsize of the tick labels
    plt.rc('legend', fontsize = MEDIUM_SIZE) # legend fontsize
    plt.rc('figure', titlesize = BIGGER_SIZE) # fontsize of the figure title
    mpl.rcParams['lines.linewidth'] = 0.5
    mpl.rcParams['lines.markersize'] = 1.5

    return figsize

def draw_figure(args, data):
    """
    Draw the figure and save it.
    """
    opt = args['doption']
    col = args['col']
    name = args['name']
    if name=='default':
        name = 'gemmes_'+opt+'.pdf'

    if opt=='all':
        nbs = 35
        nbr, nbc = nbs//5, nbs//7
        datp = data[:,1:]
        Vars = VARS[1:]
    else:
        nbr, nbc = 1, 1
        datp = data[:,col:col+1]
        Vars = VARS[col:col+1]
    time = data[:,0]

    # make figure
    figsize = change_mpl_opt(opt)
    fig, axes = plt.subplots(nrows=nbr, ncols=nbc, figsize=figsize)
    j = -1
    if len(Vars)==1:
        axes = np.asarray([[axes]])
    for n, var in enumerate(Vars):
        if n%5==0:
            j += 1
            i = 0
        ax = axes[j, i]
        try:
            ax.plot(time, datp[:,n])
            ax.set_xlabel('t')
            ax.set_ylabel(var)
        except IndexError:
            continue
        i += 1

    if args['plot']:
        plt.show()
    else:
        print('\n> Writing file {}.'.format(name))
        plt.tight_layout()
        plt.savefig(fname=name)
    plt.close(fig)

def ferror(key=0, msg=''):
    """
    Error function.
    """
    if key==0:
        print('\n> Post process run with sucess.\n')
    else:
        print('\n> Error {:d}'.format(key))
        if key==1:
            print('>     Need arguments in option {}, try -h'.format(msg))

        elif key==2:
            print('>     doption must be in the list, try -h.')

        print('\n> Post process run with an error.\n')

    # exit script
    print()
    sys.exit()

# global constants ------------------------------------------------------------
VARS = ['all', 'capital', 'npop', 'debt', 'wage', 'productivity', 'price',
        'eland', 'sigma', 'gsigma', 'co2at', 'co2up', 'co2lo', 'temp', 
        'temp0', 'pbs', 'pcar', 'omega', 'lambda', 'debtratio', 'gdp0', 'gdp',
        'eind', 'inflation', 'abat', 'n_red_fac', 'smallpi', 'smallpi_k',
        'dam', 'dam_k', 'dam_y', 'fexo', 'find', 'rcb']
ARGD = {'name':'default', 'doption':'all', 'col':0, 'plot':False}
FIGSIZE = (20, 20)
SIZE = 8

# script ----------------------------------------------------------------------
if __name__=='__main__':
    print('\n> Script is launched')
    args = load_args()
    data = load_data()
    draw_figure(args, data)
    ferror()
