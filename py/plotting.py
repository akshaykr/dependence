import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, mlab
from matplotlib import font_manager


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

pltcolors = ['b', 'r', 'g', 'k', 'm']
styles = ["-", "--", ":"]
markerstyles = ['s', 'o']

def plot_log_log(s, d=1, savefig=False):
    names = ["plugin", "linear"]
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ls = []
    ls2 = []
    line_names = []
    for i in range(len(names)):
        try:
            name = names[i]
            f = open("./data/%s_error_d=%d_s=%s.out" % (name, d, s)).readlines()
            ns = [float(x) for x in f[0].split(" ")[1:] if float(x) < 2000]
            ms = [float(x) for x in f[1].split(" ")[1:]]
            vs = [float(x) for x in f[2].split(" ")[1:]]
            ms = ms[0:len(ns)]
            vs = vs[0:len(ns)]
#             ls.append(ax1.loglog(ns, ms, linewidth=3))
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            (line, cap_lines, barlines) = ax1.errorbar(ns, ms, vs, linewidth=2, ecolor='black', elinewidth=1)
            ls.append(line)
            line_names.append("%s" % name)
            (m,b) = np.polyfit(np.log10(ns), np.log10(ms), 1)
            print (m,b)
            ls2.append(ax1.plot(np.arange(ns[0]-1, ns[-1]+1, 10), [10**(m*np.log10(x) + b) for x in np.arange(ns[0]-1, ns[-1]+1, 10)], "--", linewidth=2))
#             line_names.append("%s fit" % name)
        except IOError:
            continue
    ax1.set_xlabel("log(n)", fontsize=16)
    ax1.set_ylabel("log(error)", fontsize=16)
    ax1.set_xlim((5, 3000))
    ax1.set_ylim((5* 10**-4, 0.5))
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax1.legend(ls, line_names, prop=font_manager.FontProperties(size=16), loc=3)
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.gcf().subplots_adjust(left=0.20)
    if savefig:
        plt.savefig("./figs/log_log_d=%d_s=%s.eps" % (d, s), type="eps", dpi=1000)

    plt.show()

def quadratic_plot(ss, savefig=False):
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ls = []
    ls2 = []
    line_names = []
    for s in ss:
        try:
            f = open("./data/quadratic_error_d=1_s=%s.out" % (s)).readlines()
            ns = [float(x) for x in f[0].split(" ")[1:] if float(x) < 2000]
            ms = [float(x) for x in f[1].split(" ")[1:]]
            vs = [float(x) for x in f[2].split(" ")[1:]]
            ms = ms[0:len(ns)]
            vs = vs[0:len(vs)]
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            (line, cap_lines, barlines) = ax1.errorbar(ns, ms, vs, linewidth=2, ecolor='black', elinewidth=1)
            ls.append(line)
            line_names.append("s = %s" % s)
            (m,b) = np.polyfit(np.log10(ns), np.log10(ms), 1)
            print m, b
#             ls2.append(ax1.plot(np.arange(np.log10(ns)[0]-1, np.log10(ns)[-1]+1), [m*x + b for x in np.arange(np.log10(ns)[0]-1, np.log10(ns)[-1]+1)], '--', linewidth=2))
            ls2.append(ax1.plot(np.arange(ns[0]-1, ns[-1]+1, 10.0), [10**(m*np.log10(x)+b) for x in np.arange(ns[0]-1, ns[-1]+1, 10.0)], '--', linewidth=2))
#             line_names.append("y = %0.2f x + c" % (m))
        except IOError:
            continue
    ax1.set_xlabel("log(n)", fontsize=14)
    ax1.set_ylabel("log(error)", fontsize=14)
    ax1.set_xlim((5, 3000))
    ax1.set_ylim((5*10**-4, 0.5))
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax1.legend(ls, line_names, prop=font_manager.FontProperties(size=14), loc=3)
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.gcf().subplots_adjust(left=0.20)
    if savefig:
        plt.savefig("./figs/quadratic_log_log_d=1.eps", type="eps", dpi=1000)

    plt.show()
    

def divergence_plot(type, ss, savefig=False):
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ls = []
    ls2 = []
    line_names = []
    for s in ss:
        try:
            if type == "l2":
                f = open("./data/l2_error_d=1_s=%s.out" % (s)).readlines()
            else:
                f = open("./data/linear_%s_error_d=1_s=%s.out" % (type, s)).readlines()
            ns = [float(x) for x in f[0].split(" ")[1:] if float(x) < 2000]
            ms = [float(x) for x in f[1].split(" ")[1:]]
            vs = [float(x) for x in f[2].split(" ")[1:]]
            ms = ms[0:len(ns)]
            vs = vs[0:len(ns)]
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            (line, cap_lines, barlines) = ax1.errorbar(ns, ms, vs, linewidth=2, ecolor='black', elinewidth=1)
            ls.append(line)
            line_names.append("s = %s" % s)
            (m,b) = np.polyfit(np.log10(ns), np.log10(ms), 1)
            print m, b
#             ls2.append(ax1.plot(np.arange(np.log10(ns)[0]-1, np.log10(ns)[-1]+1), [m*x + b for x in np.arange(np.log10(ns)[0]-1, np.log10(ns)[-1]+1)], '--', linewidth=2))
            ls2.append(ax1.plot(np.arange(ns[0]-1, ns[-1]+1, 10.0), [10**(m*np.log10(x)+b) for x in np.arange(ns[0]-1, ns[-1]+1, 10.0)], '--', linewidth=2))
#             line_names.append("y = %0.2f x + c" % (m))
        except IOError:
            continue
    ax1.set_xlabel("log(n)", fontsize=14)
    ax1.set_ylabel("log(error)", fontsize=14)
    ax1.set_xlim((5, 3000))
    ax1.set_ylim((5*10**-4, 0.5))
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax1.legend(ls, line_names, prop=font_manager.FontProperties(size=14), loc=3)
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.gcf().subplots_adjust(left=0.20)
    if savefig:
        plt.savefig("./figs/%s_log_log_d=1.eps" % (type), type="eps", dpi=1000)

    plt.show()
        
