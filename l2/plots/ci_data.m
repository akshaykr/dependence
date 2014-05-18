function [ns, ps] = ci_data(ns, d, iters)

ps = [];
theta = 1/(2*sqrt(pi))^d * 2*(exp(1/4) -1)/exp(1/4);
for n=ns,
    success = 0;
    for i=1:iters,
        if d == 1,
            x = normrnd(0, 1, [d,n]);
            y = normrnd(1, 1, [d,n]);
        else,
            v = repmat(0, 1, d);
            x = mvnrnd(v, eye(d), n)';
            v(1) = 1;
            y = mvnrnd(v, eye(d), n)';
        end;
        [l2 lb ub] = confidence_interval(x,y,0.05);
        if theta >= lb && theta <= ub,
            success = success + 1;
        end;
    end;
    ps = [ps success/iters];
    fprintf('n = %d, p = %0.3f\n', n, success/iters);
end;
