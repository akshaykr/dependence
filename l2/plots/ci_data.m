function [ns, ps] = ci_data(ns, d, iters, split, ci_method)
%% evaluate the quality of a confidence interval for the L_2 divergence estimator
%% in d dimensions where data is drawn from two gaussians. 
%% split is a boolean that specifies whether we should be data-splitting on the
%% estimator or not. Usually you should set split = 1 and ci_method = 1. 
%% Returns a vector of samples sizes (which you passed in) and the empirical
%% probabilities of trapping the true parameter (averaged over the number of
%% iterations).

ps = [];
theta = 2/(2*sqrt(pi))^d * (1 - exp(-d/4));
for n=ns,
    success = 0;
    for i=1:iters,
        if d == 1,
            x = normrnd(0, 1, [d,n]);
            y = normrnd(1, 1, [d,n]);
        else,
            v = repmat(0, 1, d);
            x = mvnrnd(v, eye(d), n)';
            v = repmat(1, 1, d);
            y = mvnrnd(v, eye(d), n)';
        end;
        [l2 lb ub] = confidence_interval(x,y,0.10, split, ci_method);
        %% fprintf('truth %0.2e l2 %0.2e lb %0.2e ub %0.2e\n', theta, ...
        %%       l2, lb, ub);
        if theta >= lb && theta <= ub,
            success = success + 1;
        end;
    end;
    ps = [ps success/iters];
    fprintf('n = %d, p = %0.3f\n', n, success/iters);
end;
