function [l2] = kernel_l2(X,Y,varargin)
   p = inputParser;
   p.addParamValue('beta', 2, @isscalar);
   p.parse(varargin{:});
   Prms = p.Results;

   n1 = size(X,2);
   n2 = size(Y,2);
   d = size(X,1);

   h = min(n1,n2)^(-2/(4*Prms.beta+d));

   T1 = GaussKernel(h, X');
   T2 = GaussKernel(h, Y');
   T3 = GaussKernel(h, X', Y');
   
   T1 = 1/(n1*(n1-1)) * sum(sum(T1 - diag(diag(T1))));
   T2 = 1/(n2*(n2-1)) * sum(sum(T2 - diag(diag(T2))));
   T3 = 2/(n1*n2) * sum(sum(T3));
   l2 = T1 + T2 - T3;
   