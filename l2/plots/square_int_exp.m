%% Script for looking at the error and the confidence interval of 
%% the square integral estimator.
addpath('../common/');

d = 3;
ns = 100:100:1000;
iters = 100;

theta = (1/(2*sqrt(pi)))^d;
asymp_var = 4*( (2*sqrt(3)*pi)^(-d) - (2*sqrt(pi))^(-2*d));

zalpha = norminv(0.1/2, 0, 1);
errs = [];
ps = [];

for n=ns,
    thats = zeros(iters, 1);
    success = 0;
    for i=1:iters,
        X = normrnd(0, 1, [d,n]);
        that = square_int_estimator(X);
        thats(i) = abs(that-theta);
        lb = that + sqrt(asymp_var)*zalpha/sqrt(n);
        ub = that - sqrt(asymp_var)*zalpha/sqrt(n);
        %% fprintf('true=%0.3f est=%0.3f lb=%0.3f ub=%0.3f\n', theta, ...
        %%         that, lb, ub);
        if (theta >= lb && theta <= ub),
            success = success + 1;
        end;
    end;
    fprintf('n=%d, err=%0.2e, p=%0.2f\n', n, mean(thats), success/iters);
    errs = [errs mean(thats)];
    ps = [ps success/iters];
end;
