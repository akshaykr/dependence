function [theta, vals] = bootstrap(X,Y)
% Perform the bootstrap on the projection estimator to come up with
% the estimator and the empirical bootstrap distribution.

  n1 = size(X,2);
  n2 = size(Y,2);
  d = size(X,1);
  theta = projDivergence(X,Y);
  B = 1000;
  [bsx, bootsamx] = bootstrp(B, [], X);
  [bsy, bootsamy] = bootstrp(B, [], Y);

  vals = [];
  for i=1:B,
      Xsub = X(1:d, bootsamx(1:n1,i));
      Ysub = Y(1:d, bootsamy(1:n2,i));
      vals(i) = projDivergence(Xsub,Ysub);
  end;