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
    """
    Plot the error of the estimators on a log-log scale.
    We attempt to plot errors of "plugin" and "linear" estimators.
    s -- the smoothness we want to plot.
    """
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
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            (line, cap_lines, barlines) = ax1.errorbar(ns, ms, vs, linewidth=2, ecolor='black', elinewidth=1)
            ls.append(line)
            line_names.append("%s" % name)
            (m,b) = np.polyfit(np.log10(ns), np.log10(ms), 1)
            print (m,b)
            ls2.append(ax1.plot(np.arange(ns[0]-1, ns[-1]+1, 10), [10**(m*np.log10(x) + b) for x in np.arange(ns[0]-1, ns[-1]+1, 10)], "--", linewidth=2))
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
    """
    Plot the rate of convergence of the quadratic estimator for several smoothness values.
    """
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
            ls2.append(ax1.plot(np.arange(ns[0]-1, ns[-1]+1, 10.0), [10**(m*np.log10(x)+b) for x in np.arange(ns[0]-1, ns[-1]+1, 10.0)], '--', linewidth=2))
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
    """
    Plot the rates of convergence for a specific divergence for several smoothness values.
    """
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
            ls2.append(ax1.plot(np.arange(ns[0]-1, ns[-1]+1, 10.0), [10**(m*np.log10(x)+b) for x in np.arange(ns[0]-1, ns[-1]+1, 10.0)], '--', linewidth=2))
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
        
def plot_estimator_rate(ns, ms, vs, gamma):
    """
    Generate three plots. 
    One is an error bar plot of the error of the estimator.
    The next is n vs \exp\{\log(err(n)) + \gamma \log n\} which should pull out the constant in the rate.
    The third is n vs -\log(err(n)/\log(n) which should approach \gamma
    """
    fig = plt.figure(figsize=(15, 5))
    ax1 = fig.add_subplot(131)
    ax1.errorbar(ns, ms, vs)
    ax1.set_xlabel("Number of samples (n)")
    ax1.set_ylabel("Error |That - T|")
    ax1.set_xscale('log')
    ax2 = fig.add_subplot(132)
    ax2.plot(ns, [ms[i]*ns[i]**gamma for i in range(len(ns))])
    ax2.set_xlabel("Number of samples (n)")
    ax2.set_ylabel("Error*n^{\gamma}")
    ax3 = fig.add_subplot(133)
    ax3.plot(ns, [-np.log(ms[i])/np.log(ns[i]) for i in range(len(ns))])
    ax3.set_xlabel("Number of samples (n)")
    ax3.set_ylabel("-log(error)/log(n)")
    plt.show()

def plot_data_log_log(ns, ms, fig=None):
    """
    Plot the data on a log-log plot.
    This differs from plot_log_log because it takes in ns and ms rather than the files.
    """
    if fig == None:
        fig = plt.figure(figsize=(5,5))
    ax1 = fig.add_subplot(111)
    ax1.plot(np.log(ns), np.log(ms))
    [m,b] = np.polyfit(np.log(ns), np.log(ms), 1)
    ax1.plot(np.arange(np.log(ns)[0]-1, np.log(ns)[-1]+1), [m*x + b for x in np.arange(np.log(ns)[0]-1, np.log(ns)[-1]+1)])
    ax1.set_xlabel("log(n)")
    ax1.set_ylabel("log(error)")
    plt.show()
    return (m,b)

def plot_from_file(file, gamma, fig=None):
    """
    Plot the data from file. Generate the three plots of the rates of convergence.
    """
    f = open(file).readlines()
    ns = [float(x) for x in f[0].split(" ")[1:]]
    ms = [float(x) for x in f[1].split(" ")[1:]]
    vs = [float(x) for x in f[2].split(" ")[1:]]

    plot_estimator_rate(ns, ms, vs, gamma)
    (m,b) = plot_log_log(ns, ms, fig=fig)
    return (m,b)

def compare(s):
    """
    Compare the linear and plugin estimators for a set of ss.
    We have another routine for this now.
    """
    est_types = ["plugin", "linear"]
    ns = []
    ms = []
    vs = []
    for est_type in est_types:
        f = open("./data/%s_error_d=1_s=%s.out" % (est_type, s)).readlines()
        ns.append([float(x) for x in f[0].split(" ")[1:]])
        ms.append([float(x) for x in f[1].split(" ")[1:]])
        vs.append([float(x) for x in f[2].split(" ")[1:]])
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    lines = [ax.plot(ns[i], ms[i]) for i in range(len(est_types))]
    ax.legend(lines, ["%s" % (est_type) for est_type in est_types], prop=font_manager.FontProperties(size=14), loc=1)
    ax.set_xlabel("n")
    ax.set_ylabel("Error")
    ax.set_xscale("log")


