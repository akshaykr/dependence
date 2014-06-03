function [ns, ps] = ci_data(ns, d, iters, split, ci_method)

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
