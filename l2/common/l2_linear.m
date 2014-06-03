function [l2] = l2_linear(X,Y,varargin)
% Linear time estimator of the L2 divergence.
% Splits the data into two groups (2/3, 1/3) and on the first group 
% Estimates \int p^2(x) with \sum_{i=1}^{n/3} K(X_{2i-1}, X_{2i})
% Estimates \int q^2(x) with \sum_{i=1}^{n/3} K(Y_{2i-1}, Y_{2i})
% Estimates \int pq with \sum_{i=1}^{n/3} K(X_i, Y_i)

   p = inputParser;
   p.addParamValue('beta', 4, @isscalar);
   p.addParamValue('debug',0, @isscalar);
   p.addParamValue('fast', 0, @isscalar);
   p.parse(varargin{:});
   Prms = p.Results;

   n1 = size(X,2);
   n2 = size(Y,2);
   d = size(X,1);

   X1 = X(:,1:2*n1/3);
   X2 = X(:,2*n1/3:end);
   Y1 = Y(:,1:2*n2/3);
   Y2 = Y(:,2*n1/3:end);

   