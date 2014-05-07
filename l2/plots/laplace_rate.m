function [ns, ms, vs] = laplace_rate(ns, mu, s, iters),
%% numerically integrate the difference between the two laplace
%% distributions.

%% P_1 is L(-mu/2, 1) while P_2 is L(mu/2, s).

pts = -10:0.01:10;
theta = mean((1/2*exp(- abs(pts + mu/2)) - 1/(2*s)*exp(-abs(pts - ...
                                                  mu/2)/(s))).^2);
fprintf('theta = %0.4f\n', theta);

ms = [];
vs = [];
for n=ns,
    scores = [];
    for i=1:iters,
        xpre = unifrnd(-1/2, 1/2, [1,n]);
        x = -mu/2 - sign(xpre) .* log(1 - 2*abs(xpre));
        ypre = unifrnd(-1/2, 1/2, [1,n]);
        y = mu/2 - sign(ypre) .* log(1 - 2*abs(ypre));
        [l2] = kernel_l2(x,y);
        scores = [scores l2];
    end;
    ms = [ms mean(scores)];
    vs = [vs var(scores)];
    fprintf('n = %d, m = %0.3f, v = %0.3f\n', n, mean(scores), var(scores));
end;
return;