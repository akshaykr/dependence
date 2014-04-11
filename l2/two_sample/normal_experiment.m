function [psuccess] = normal_experiment(mus, test),
%% look at the probability of rejecting the null hypothesis when
%% data is drawn from n(0,1) and n(1,1) across the possible samples.
%% use the test statistic @test.

psuccess = [];
n = 100;

for mu=mus,
    s = 0;
    for i=1:100,
        x = normrnd(0, 1, [1,n]);
        y = normrnd(mu, 1, [1,n]);
        [h,p,that] = test(x,y);
        s = s + h;
    end;
    fprintf('mu=%0.2f, s = %d\n', mu, s);
    psuccess = [psuccess s];
end;