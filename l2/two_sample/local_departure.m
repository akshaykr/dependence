function [psuccess] = local_departure(mus, test)

psuccess = [];
iters = 100;
n = 20;

for mu=mus,
    s = 0;
    for i=1:iters,
        x = normrnd(0, 2, [1,n]);

        %% mixture of normals
        z = binornd(1, 0.5, [1,n]);
        y = z.*normrnd(mu, 1, [1,n]) + (1-z).*normrnd(-mu, 1, [1,n]);
        %% different variance normal
        % y = normrnd(0, mu, [1,n]);

        %% laplace different variance. 
        % ypre = unifrnd(-1/2, 1/2, [1,n]);
        % y = -mu*sign(ypre) .* log(1 - 2* abs(ypre));
        %% laplace different mean.
        % y = mu - sign(ypre) .* log(1-2*abs(ypre));

        %% cauchy different variance.
        % ypre = unifrnd(-pi/2, pi/2, [1,n]);
        % y = tan(ypre)*mu;
        %% cauchy different mean.
        % y = mu + tan(ypre);
        [h,p,that] =  test(x,y);
        s = s+h;
    end;
    fprintf('mu=%0.2f s = %d\n', mu, s);
    psuccess = [psuccess s/iters];
end;