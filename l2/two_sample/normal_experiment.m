function [psuccess] = normal_experiment(mus, test, d),
%% look at the probability of rejecting the null hypothesis when
%% data is drawn from n(0,1) and n(mu,1) across the possible samples.
%% use the test statistic @test.

psuccess = [];
n = 40;
iters = 200;

for mu=mus,
    s = 0;
    for i=1:iters,
        if d == 1,
            x = normrnd(0, 1, [d,n]);
            y = normrnd(1, 1, [d,n]);
        else,
            v = repmat(0, 1, d);
            x = mvnrnd(v, eye(d), n)';
            v = repmat(1, 1, d);
            y = mvnrnd(mu*v, eye(d), n)';
        end;
        [h,p,that] = test(x,y);
        s = s + h;
    end;
    fprintf('mu=%0.2f, s = %d\n', mu, s);
    psuccess = [psuccess s/iters];
end;
