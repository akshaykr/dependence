addpath('../common/');

ns = 50:50:1000;
iters = 50;
d = 1;


if d == 1
    l2_dif = @(x) (normpdf(x, 0.5, 1) - normpdf(x, -0.5, 1)).^2;
    theta = quad(l2_dif, -1000, 1000);
else,
    l2_dif = @(x) (mvnpdf(x,repmat(-0.5,1, d),eye(d)) - mvnpdf(x,repmat(0.5,1, ...
                                                      d),eye(d))).^2;
    %% Perform Monte Carlo integration.
    pts = 10000000;
    draws = unifrnd(-2, 2, [pts d]);
    theta = mean(l2_dif(draws));
end;

% theta = 0;
fprintf('theta = %0.3f\n', theta);

ms = [];
vs = [];

for n=ns,
    scores = [];
    for i=1:iters,
        if d == 1,
            x = normrnd(-0.5, 1, [1 n]);
            y = normrnd(0.5, 1, [1 n]);
        else,
            x = mvnrnd(repmat(-0.5,1,d), eye(d), n)';
            y = mvnrnd(repmat(0.5,1,d), eye(d), n)';
        end;
        [l2] = kernel_l2(x,y);
        scores = [scores abs(l2 - theta)];
    end;
    ms = [ms mean(scores)];
    vs = [vs var(scores)];
    fprintf('n = %d, m = %0.3f, v = %0.3f\n', n, mean(scores), var(scores));
end;
