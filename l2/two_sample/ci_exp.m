addpath('../common/');

ns = 100:100:800;
boot = 0;

l2_dif = @(x) (normpdf(x,0,1) - normpdf(x,1,1)).^2;

theta = quad(l2_dif, -1000, 1000);
fprintf('theta = %0.3f\n', theta);

iters = 50;
B = 100;
pvar = [];
pboot = [];
bias = [];
for n=ns,
    success = 0;
    b_success = 0;
    biases = [];
    for i=1:iters,
        x = normrnd(0, 1, [1,n]);
        y = normrnd(1, 1, [1,n]);

        [l2 lb ub] = confidence_interval(x,y,0.05);
        biases = [biases l2 - theta];
        if theta >= lb && theta <= ub,
            success = success + 1;
        end;
        if boot == 1,
            [l2 vals lb ub] = bootstrap(x,y,0.05, B);
            if theta >= lb && theta <= ub,
                b_success = b_success + 1;
            end;
        end;
    end;
    bias = [bias mean(biases)];
    pvar = [pvar success/iters];
    pboot = [pboot -1];
    fprintf('n = %d, bias = %0.3f p_var = %0.3f p_boot = %0.3f\n', n, mean(biases), success/iters, ...
            b_success/iters);
end;
