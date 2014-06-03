%% Basic synthetic spectral clustering example on a hierarchical model.
%% Cluster marginals are gaussians with mean (0, \ldots, 0) and \mu (1, \ldots, 1)
%% both in d-dimensions and variance sigma. We draw n distributions per cluster marginal
%% m samples per distribution. Then we compute the gram matrix of L_2 divergences between these
%% distributions and run (normalized) spectral clustering to obtain a solution. 
%% We repeat this for different values of m and mu and output a matrix of error rates across both dimensions.

addpath('../common/');

n = 50;
ms = 10:10:200;
mus = 0.25:0.25:1.0;
d = 5;
sigma = 0.1;
iters = 10;

errs = zeros(size(mus, 2), size(ms,2));
for j=1:size(mus,2),
    for k = 1:size(ms,2),
        err = [];
        for i=1:iters,
            [Data, Labels] = mixture_model(n,ms(k),d,mus(j),sigma);
            [c1, c2] = cluster_mixture(Data);
            s1 = size(setdiff(c1, find(Labels == 1)),1);
            s2 = size(setdiff(c1, find(Labels == -1)),1);
            err = [err min(s1, s2)/n];
        end
        fprintf('mu = %0.2f m = %d err = %0.2f\n', mus(j), ms(k), mean(err));
        errs(j,k) = mean(err);
    end;
end;
