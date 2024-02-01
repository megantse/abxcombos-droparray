# import packages
import matplotlib.pyplot as plt
from cycler import cycler

# define plotting parameters
def parameters():
    fontsize = 14
    plt.rcParams['axes.spines.right']=False
    plt.rcParams['axes.spines.top']=False
    plt.rcParams['axes.linewidth']=3
    plt.rcParams['axes.labelsize']=fontsize
    plt.rcParams['lines.linewidth']=2
    plt.rcParams['xtick.labelsize']=fontsize
    plt.rcParams['ytick.labelsize']=fontsize
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['font.size']=fontsize
    plt.rcParams['xtick.major.width']=1.5
    plt.rcParams['ytick.major.width']=1.5
    plt.rcParams['contour.negative_linestyle'] = 'solid'

    plt.rcParams['savefig.bbox']='Tight'
    plt.rcParams['savefig.dpi']=300
    #
    CB_color_cycle = ["#b65a78",
                    "#27b5a2",
                    "#e28776",
                    "#00515e",
                    "#fff6b9",
                    "#543761",
                    "#afa156",
                    "#ffeef8",
                    "#7c772e"]

    plt.rcParams['axes.prop_cycle'] = cycler('color', CB_color_cycle)
