from helper import fast_integration
import kde
import numpy as np

"""
A small experiment to check if the 1-density version of our linear estimator is asymptotically normal.
"""

def test_uni_clt(P, alpha, n1, n2, s, iters=200):
    """
    Repeatedly sample data from the density and compute the 1-density version of the linear estimator.
    Return a vector of length iters, with the error for each trial
    """
    vals = []
    truth = fast_integration(lambda x: np.power(P.eval(np.matrix(x)), alpha), [0.0], [1.0])
    for i in range(iters):
        d1 = P.sample(n1)
        d2 = P.sample(n2)
        k = kde.KDE(d1, s)
        score = np.mean(np.power(k.eval(d2), alpha-1))
        vals.append(score - truth)
#         print "i = %d, score = %0.2f" % (i, score-truth)
    return np.array(vals)


if __name__=="__main__":
    import density
    import matplotlib.pyplot as plt
    P = density.UniTrigDensity(2, 1)
    scores = []
    vals = range(50, 401, 50)
    for n in vals:
        print "%d samples" % (n)
        scores.append(test_uni_clt(P, 0.5, n/2, n/2, 2))

    fig = plt.figure(figsize=(15,5))
    ax1 = fig.add_subplot(141)
    [ax1.hist(scores[i], bins=50) for i in range(0, len(scores), 2)]
    ax1.set_title("Un-normalized histograms")

    ax1 = fig.add_subplot(142)
    [ax1.hist(np.sqrt(vals[i])*scores[i], bins=50) for i in range(0, len(scores), 2)]
    ax1.set_title("Normalized histograms")

    ax1 = fig.add_subplot(1,4,3)
    ax1.plot(vals, [np.var(scores[i]) for i in range(len(scores))])
    ax1.set_title("Variance vs samples")

    ax1 = fig.add_subplot(1,4,4)
    ax1.plot(vals, [np.var(np.sqrt(vals[i])*scores[i]) for i in range(len(scores))])
    ax1.set_title("Normalized variance vs Samples")
    plt.show()
