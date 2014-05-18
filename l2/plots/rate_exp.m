addpath('../common/');

ns = 20:20:1000;
iters = 200;

d = 1;
%% OLD -- I'm looking at the p=q situation so theta = 0.
if d == 1
    l2_dif = @(x) (normpdf(x, 0.5, 1) - normpdf(x, -0.5, 1)).^2;
    theta = quad(l2_dif, -1000, 1000);
else,
    l2_dif = @(x) (mvnpdf(x,repmat(-0.5,1, d),eye(d)) - mvnpdf(x,repmat(0.5,1, ...
                                                      d),eye(d))).^2;
    %% Perform Monte Carlo integration.
    pts = 10000;
    draws = unifrnd(-2, 2, [pts d]);
    theta = mean(l2_dif(draws));
end;

theta = 0;
fprintf('theta = %0.3f\n', theta);
fprintf('d = 1\n');
[ns, ms1, vs1] = rate_data(ns, 1, iters);
fprintf('d = 2\n');
[ns, ms2, vs2] = rate_data(ns, 2, iters);
fprintf('d = 3\n');
[ns, ms3, vs3] = rate_data(ns, 3, iters);

save('./rate_data.mat', 'ns', 'ms1', 'vs1', 'ms2', 'vs2', 'ms3', 'vs3');
