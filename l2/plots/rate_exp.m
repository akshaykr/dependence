addpath('../common/');

ns = 50:50:1000;
iters = 50;
d = 10;

l2_dif = @(x) (normpdf(x,0,1) - normpdf(x,1,1)).^2;
theta = quad(l2_dif, -1000, 1000);

% theta = 0;
fprintf('theta = %0.3f\n', theta);

ms = [];
vs = [];

for n=ns,
    scores = [];
    for i=1:iters,
        x = normrnd(0, 1, [1,n]);
        y = normrnd(1, 1, [1,n]);

        [l2] = kernel_l2(x,y);
        scores = [scores abs(l2 - theta)];
    end;
    ms = [ms mean(scores)];
    vs = [vs var(scores)];
    fprintf('n = %d, m = %0.3f, v = %0.3f\n', n, mean(scores), var(scores));
end;
