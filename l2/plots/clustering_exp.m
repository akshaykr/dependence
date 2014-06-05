%% Basic synthetic spectral clustering example on a hierarchical model.
%% Cluster marginals are gaussians with mean (0, \ldots, 0) and \mu (1, \ldots, 1)
%% both in d-dimensions and variance sigma. We draw n distributions per cluster marginal
%% m samples per distribution. Then we compute the gram matrix of L_2 divergences between these
%% distributions and run (normalized) spectral clustering to obtain a solution. 
%% We repeat this for different values of m and mu and output a matrix of error rates across both dimensions.

addpath('../common/');
addpath('../../npdiv_export/');
addpath('../../npdiv_export/lib/');
n = 20;
ms = 10:10:100;
mus = 0.25:0.25:1.0;
d = 3;
sigma = 0;
iters = 2;

kde_errs = zeros(size(mus, 2), size(ms,2));
nn_errs = zeros(size(mus, 2), size(ms,2));
for j=1:size(mus,2),
    for k = 1:size(ms,2),
        kde_err = [];
        nn_err = [];
        for i=1:iters,
            [Data, Labels] = mixture_model(n,ms(k),d,mus(j),sigma);
            [c1, c2] = cluster_mixture(Data, 0);
            s1 = size(setdiff(c1, find(Labels == 1)),1) + size(setdiff(c2, find(Labels == -1)), 1);
            s2 = size(setdiff(c1, find(Labels == -1)),1) + size(setdiff(c2, find(Labels == 1)), 1);
            kde_err = [kde_err min(s1, s2)/(2*n)];
            [c1, c2] = cluster_mixture(Data, 1);
            s1 = size(setdiff(c1, find(Labels == 1)),1) + size(setdiff(c2, find(Labels == -1)), 1);
            s2 = size(setdiff(c1, find(Labels == -1)),1) + size(setdiff(c2, find(Labels == 1)), 1);
            nn_err = [nn_err min(s1, s2)/(2*n)];
        end
        fprintf('mu = %0.2f m = %d kde_err = %0.2f nn_err = %0.2f\n', mus(j), ms(k), mean(kde_err), mean(nn_err));
        kde_errs(j,k) = mean(kde_err);
        nn_errs(j,k) = mean(nn_err);
    end;
end;
