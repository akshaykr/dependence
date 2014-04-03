addpath('../common/');

ns = 10:10:200;

l2_dif = @(x) (normpdf(x,0,1) - normpdf(x,1,1)).^2;

theta = quad(l2_dif, -1000, 1000);
fprintf('theta = %0.3f\n', theta);

iters = 100;
B = 100;
pvar = [];
pboot = [];
for n=ns,
    success = 0;
    b_success = 0;
    for i=1:iters,
        x = normrnd(0, 1, [1,n]);
        y = normrnd(1, 1, [1,n]);

        [l2 lb ub] = confidence_interval(x,y,0.05);
        if theta >= lb && theta <= ub,
            success = success + 1;
        end;
        [l2 vals lb ub] = bootstrap(x,y,0.05, B);
        if theta >= lb && theta <= ub,
            b_success = b_success + 1;
        end;
    end;
    pvar = [pvar success/iters];
    pboot = [pboot b_success/iters];
    fprintf('n = %d, p_var = %0.3f p_boot = %0.3f\n', n, success/iters, ...
            b_success/iters);
end;
