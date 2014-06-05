function [Data, Labels] = mixture_model(n,m,d,mu,sigma),
%% Two level mixture model in d dimensions.
%% class marginals are N(0, I), N(\mu, I) and this is the prior
%% distribution on the individual distributions means. 
%% draw n points per class and m samples per point.

Data = {};
% Data = zeros(m, d, 2*n);
Labels = zeros(2*n, 1);

if sigma == 0,
    c0_means = zeros(n,d);
else,
    c0_means = mvnrnd(zeros(1,d), sigma*eye(d), n);
end
for i=1:n,
    x = mvnrnd(c0_means(i,:), eye(d), m);
    Data{i} = x;
    Labels(i) = -1;
end;

if sigma == 0,
    c1_means = mu*ones(n,d);
else
    c1_means = mvnrnd(mu*ones(1,d), sigma*eye(d), n);
end;
for i=1:n,
    x = mvnrnd(c1_means(i,:), eye(d), m);
    Data{i+n} = x;
    Labels(i+n) = 1;
end;
