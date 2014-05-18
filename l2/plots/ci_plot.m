addpath('../common/');

n = 50;

l2_dif = @(x) (normpdf(x,0,1) - normpdf(x,1,1)).^2;

theta = quad(l2_dif, -100, 100);
fprintf('theta = %0.2f\n', theta);

iters = 100;
l2s = [];
lbs = [];
ubs = [];
success = 0;
for i=1:iters,
    x = normrnd(0,1,[1, n]);
    y = normrnd(1,1,[1, n]);


    [l2, lb,ub] = confidence_interval(x,y,0.05);
    l2s = [l2s l2];
    lbs = [lbs lb];
    ubs = [ubs ub];
    fprintf('iter = %d, l2 = %0.2f, lb = %0.2f, ub = %0.2f\n', i, l2, lb, ub);
    if theta >= lb && theta <= ub,
        success = success + 1;
    end;
end;

f = figure();
plot(1:1:iters, repmat(theta, [1 iters]));
hold on;
for i=1:iters,
    plot([i i], [lbs(i), ubs(i)]);
end;

pest = success/iters;