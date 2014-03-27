function [theta, vals] = bootstrap_square(X, varargin)
% Perform the bootstrap on the projection estimator to come up with
% the estimator and the empirical bootstrap distribution.
  B = 100;
  if size(varargin,1) == 1,
      B = varargin{1};
  end;

  n1 = size(X,2);
  d = size(X,1);
  theta = square_int_estimator(X);
  [bsx, bootsamx] = bootstrp(B, [], X);

  vals = [];
  for i=1:B,
      Xsub = X(1:d, bootsamx(1:n1,i));
      vals(i) = square_int_estimator(Xsub);
  end;