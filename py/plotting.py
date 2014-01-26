import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

pltcolors = ['b', 'r', 'g', 'k', 'm']
styles = ["-", "--", ":"]

def plot_log_log(s, d=1, savefig=False):
    names = ["plugin", "linear", "quadratic"]
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ls = []
    line_names = []
    for name in names:
        try:
            f = open("./data/%s_error_d=%d_s=%s.out" % (name, d, s)).readlines()
            ns = [float(x) for x in f[0].split(" ")[1:]]
            ms = [float(x) for x in f[1].split(" ")[1:]]
            vs = [float(x) for x in f[2].split(" ")[1:]]
            ls.append(ax1.plot(np.log(ns), np.log(ms)))
            line_names.append("%s" % name)
            (m,b) = np.polyfit(np.log(ns), np.log(ms), 1)
            print (m,b)
            ls.append(ax1.plot(np.arange(np.log(ns)[0]-1, np.log(ns)[-1]+1), [m*x + b for x in np.arange(np.log(ns)[0]-1, np.log(ns)[-1]+1)], "--"))
            line_names.append("%s fit" % name)
        except IOError:
            continue
    ax1.set_xlabel("log(n)")
    ax1.set_ylabel("log(error)")
    ax1.legend(ls, line_names, prop=font_manager.FontProperties(size=14), loc=1)
    plt.gcf().subplots_adjust(bottom=0.20)
    if savefig:
        plt.savefig("./figs/log_log_d=%d_s=%s.eps" % (d, s), type="eps", dpi=1000)

    plt.show()

    
