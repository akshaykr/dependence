import matplotlib.pyplot as plt
import numpy as np
from matplotlib import (
    font_manager,
    rc,
)

"""
Plot the theoretical rates of convergence of the estimators in the paper.
"""

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

pltcolors = ['b', 'r', 'g', 'k', 'm']
styles = ["-", "--", ":"]

if __name__=="__main__":
    d = 10
    f1 = lambda s:  float(4*s)/(4*s+d)
    f2 = lambda s:  float(3*s)/(2*s+d)
    f3 = lambda s:  1.0/2
    f4 = lambda s:  float(2*s)/(2*s+d)
    f5 = lambda s:  float(1*s)/(2*s+d)

    x = np.arange(0, d, 0.01)
    fig = plt.figure(figsize=(20, 3.5))
    ax1 = fig.add_subplot(141)

    l1 = ax1.plot(x, [f1(i) for i in x], linewidth=3, linestyle=styles[2], color=pltcolors[2])
    l2 = ax1.plot(x, [f2(i) for i in x], linewidth=3, linestyle=styles[0], color=pltcolors[0])
    l3 = ax1.plot(x, [f3(i) for i in x], linewidth=3, linestyle=styles[1], color=pltcolors[1])
    plt.legend([l2, l3, l1], ["remainder", "linear", "quadratic"], prop=font_manager.FontProperties(size=14), loc=4)
    ax1.set_xlabel('Smoothness (s)', fontsize=14)
    ax1.set_ylabel('Rate of convergence ($\gamma$)')

    ax2 = fig.add_subplot(142)
    l1 = ax2.plot(x, [f4(i) for i in x], linewidth=3, linestyle=styles[0], color=pltcolors[0])
    l2 = ax2.plot(x, [f3(i) for i in x], linewidth=3, linestyle=styles[1], color=pltcolors[1])
    plt.legend([l1, l2], ["remainder", "linear"], prop=font_manager.FontProperties(size=14), loc=4)
    ax2.set_xlabel('Smoothness (s)', fontsize=14)
    ax2.set_ylabel('Rate of convergence ($\gamma$)')

    ax3 = fig.add_subplot(143)
    l1 = ax3.plot(x, [f5(i) for i in x], linewidth=3, linestyle=styles[0], color=pltcolors[0])
    plt.legend([l1], ["remainder"], prop=font_manager.FontProperties(size=14), loc=4)
    ax3.set_xlabel('Smoothness (s)', fontsize=14)
    ax3.set_ylabel('Rate of convergence ($\gamma$)')

    ax4 = fig.add_subplot(144)
    l1 = ax4.plot(x, [f1(i) for i in x], linewidth=3, linestyle=styles[2], color=pltcolors[2])
    l3 = ax4.plot(x, [f3(i) for i in x], linewidth=3, linestyle=styles[1], color=pltcolors[1])
    plt.legend([l3, l1], ["linear", "quadratic"], prop=font_manager.FontProperties(size=14), loc=4)
    ax4.set_xlabel('Smoothness (s)', fontsize=14)
    ax4.set_ylabel('Rate of convergence ($\gamma$)', fontsize=14)

    ax1.set_ylim(0, 1)
    ax2.set_ylim(0, 1)
    ax3.set_ylim(0, 1)
    ax4.set_ylim(0, 1)

    plt.gcf().subplots_adjust(bottom=0.20)

    plt.savefig("./figs/rate_plot.eps", dpi=1000, format="eps", bbox_inches="tight")

    plt.show()
