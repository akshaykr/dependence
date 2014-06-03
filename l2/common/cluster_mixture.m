function [c1 c2] = cluster_mixture(Data),

n = size(Data,3);

M = zeros(n,n);

for i=1:n,
    for j=i:n,
        M(i,j) = exp(- kernel_l2(Data(:,:,i)', Data(:,:,j)'));
        M(j,i) = M(i,j);
    end;
end;

%% project M onto PSD cone by zero-ing out negative eigenvalues.
% [u,v] = eig(M);
% v(find(v < 0)) = 0;
% M = u*v*u';

D = diag(diag(M));
% L = D - M;
Ls = D^(-1/2)*M*D^(-1/2);

[u,v] = eig(Ls);

vec = u(:,2);
c1 = find(vec > 0);
c2 = find(vec < 0);
