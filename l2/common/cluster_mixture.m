function [c1 c2] = cluster_mixture(Data),
%% Perform two way normalized spectral clustering on the set of distributions in Data
%% Data is a m x d x n tensor of m samples from each of n distributions in d dimensions. 
%% Use a similarity function based on the L_2 divergence estimator. 

n = size(Data,3);

M = zeros(n,n);

for i=1:n,
    for j=i:n,
        M(i,j) = exp(- kernel_l2(Data(:,:,i)', Data(:,:,j)'));
        M(j,i) = M(i,j);
    end;
end;


D = diag(diag(M));
Ls = D^(-1/2)*M*D^(-1/2);

[u,v] = eig(Ls);

vec = u(:,2);
c1 = find(vec > 0);
c2 = find(vec < 0);
