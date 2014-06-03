function [theta, vals, lb, ub] = bootstrap(X,Y,alpha, varargin)
% Perform the bootstrap on the kernel l2 estimator to come up with
% the estimator and the empirical bootstrap distribution.
  B = 500;
  if size(varargin,1) == 1,
      B = varargin{1};
  end;

  n1 = size(X,2);
  n2 = size(Y,2);
  d = size(X,1);
  theta = kernel_l2_no_split(X,Y);
  [bsx, bootsamx] = bootstrp(B, [], X');
  [bsy, bootsamy] = bootstrp(B, [], Y');

  vals = [];
  for i=1:B,
      Xsub = X(1:d, bootsamx(1:n1, i));
      Ysub = Y(1:d, bootsamy(1:n2, i));
      vals(i) = kernel_l2_no_split(Xsub,Ysub);
  end;

  s = sort(vals - theta);
  lb = theta - s(int64((1-alpha/2)*B));
  ub = theta - s(int64(alpha/2*B));