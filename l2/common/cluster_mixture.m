function [c1, c2, M] = cluster_mixture(Data,div),
%% Perform two way normalized spectral clustering on the set of distributions in Data
%% Data is a m x d x n tensor of m samples from each of n distributions in d dimensions. 
%% Use a similarity function based on the L_2 divergence estimator. 

% n = size(Data,3);
n = length(Data);

if div == 0,
    M = zeros(n,n);

    for i=1:n,
        for j=i:n,
            M(i,j) = exp(- kernel_l2(Data{i}', Data{j}'));
            M(j,i) = M(i,j);
        end;
    end;
else,
    M = NPDivs(Data, []);
    M = 0.5*(M+M');
    M = exp(-M);
end;


D = diag(diag(M));
Ls = D^(-1/2)*M*D^(-1/2);

[u,v] = eigs(Ls, 2);

vec = u(:,2);
c1 = find(vec > 0);
c2 = find(vec < 0);
