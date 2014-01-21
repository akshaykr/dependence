import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, mlab
from matplotlib import font_manager

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

if __name__=="__main__":
    d = 10
    f1 = lambda s:  float(4*s)/(4*s+d)
    f2 = lambda s:  float(3*s)/(2*s+d)
    f3 = lambda s:  1.0/2
    f4 = lambda s:  float(2*s)/(2*s+d)
    f5 = lambda s:  float(1*s)/(2*s+d)

    x = np.arange(0, d, 0.01)
    fig = plt.figure(figsize=(20, 4))
    ax1 = fig.add_subplot(141)

    l1 = ax1.plot(x, [f1(i) for i in x])
    l2 = ax1.plot(x, [f2(i) for i in x])
    l3 = ax1.plot(x, [f3(i) for i in x])
    plt.legend([l1, l2, l3], ["quadratic", "remainder", "linear"], prop=font_manager.FontProperties(size=14), loc=4)
    ax1.set_xlabel('Smoothness (s)', fontsize=14)
    ax1.set_ylabel('Rate of convergence ($\gamma$)')

    ax2 = fig.add_subplot(142)
    l1 = ax2.plot(x, [f4(i) for i in x])
    l2 = ax2.plot(x, [f3(i) for i in x])
    plt.legend([l1, l2, l3], ["remainder", "linear"], prop=font_manager.FontProperties(size=14), loc=4)
    ax2.set_xlabel('Smoothness (s)', fontsize=14)
    ax2.set_ylabel('Rate of convergence ($\gamma$)')

    ax3 = fig.add_subplot(143)
    l1 = ax3.plot(x, [f5(i) for i in x])
    plt.legend([l1], ["remainder"], prop=font_manager.FontProperties(size=14), loc=4)
    ax3.set_xlabel('Smoothness (s)', fontsize=14)
    ax3.set_ylabel('Rate of convergence ($\gamma$)')

    ax4 = fig.add_subplot(144)
    l1 = ax4.plot(x, [f1(i) for i in x])
    l3 = ax4.plot(x, [f3(i) for i in x])
    plt.legend([l1, l3], ["quadratic", "linear"], prop=font_manager.FontProperties(size=14), loc=4)
    ax4.set_xlabel('Smoothness (s)', fontsize=14)
    ax4.set_ylabel('Rate of convergence ($\gamma$)')

    ax1.set_ylim(0, 1)
    ax2.set_ylim(0, 1)
    ax3.set_ylim(0, 1)
    ax4.set_ylim(0, 1)

    plt.gcf().subplots_adjust(bottom=0.20)

    plt.savefig("./figs/rate_plot.eps", dpi=1000, format="eps")
    
    plt.show()
