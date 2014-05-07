function [psuccess] = laplace_experiment(mus, test),
%% look at the probability of rejecting the null hypothesis when
%% data is drawn from n(0,1) and n(mu,1) across the possible samples.
%% use the test statistic @test.
addpath('../common/');

psuccess = [];
n = 20;
iters = 100;

for mu=mus,
    s = 0;
    for i=1:iters,
        xpre = unifrnd(-1/2, 1/2, [1,n]);
        x = -mu/2 - sign(xpre) .* log(1 - 2*abs(xpre));
        ypre = unifrnd(-1/2, 1/2, [1,n]);
        y = mu/2 - sign(ypre) .* log(1 - 2*abs(ypre));
        [h,p,that] = test(x,y);
        s = s + h;
    end;
    fprintf('mu=%0.2f, s = %d\n', mu, s);
    psuccess = [psuccess s/iters];
end;
