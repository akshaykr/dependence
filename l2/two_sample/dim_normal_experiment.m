function [psuccess] = dim_normal_experiment(mu, test, ds, scaled)
%% Routine for looking at two sample testing between normals as a function
%% of dimension. 
%% Scaled indicates whether to rescale the second distribution.
%% scaled=1 makes the second distribution mu/sqrt(d) (1, \ldots, 1)
%% otherwise it is mu(1, \ldots, 1).

psuccess = [];
n = 100;
iters = 200;

for d=ds,
    s = 0;
    for i=1:iters,
        if d == 1,
            x = normrnd(0, 1, [d,n]);
            y = normrnd(mu, 1, [d,n]);
        else,
            v = repmat(0, 1, d);
            x = mvnrnd(v, eye(d), n)';
            v = repmat(1, 1, d);
            if scaled,
                y = mvnrnd(mu*v/sqrt(d), eye(d), n)';
            else,
                y = mvnrnd(mu*v, eye(d), n)';
            end;
        end;
        [h,p,that] = test(x,y);
        s = s + h;
    end;
    fprintf('d=%d, s = %d\n', d, s);
    psuccess = [psuccess s/iters];
end;
