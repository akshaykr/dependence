function [h,p,that, others] = l2_permutation(X,Y)
  that = kernel_l2_no_split(X,Y);
  n1 = size(X,2);
  n2 = size(Y,2);

  alpha = 0.05;

  B = 500;
  others = zeros(B,1);
  for b =1:B,
      inds1 = randsample(n1, n1/2);
      inds2 = randsample(n2, n2/2);
      Xnew = [X(:,inds1) Y(:,inds2)];
      Ynew = [X(:,setdiff(1:n1, inds1)) Y(:,setdiff(1:n2,inds2))];
      others(b,1) = kernel_l2_no_split(Xnew, Ynew);
  end;

  %% others = sort(others);
  %% fprintf('lb = %0.2f, ub = %0.2f, that = %0.2f\n', ...
  %%       others(floor(alpha/2*B),1), others(ceil((1-alpha/2)*B), ...
  %%                                          1), that);
  p = mean(others >= that);
  %% fprintf('that=%0.3f p=%0.3f\n', that, p);
  if p < alpha,
      h = 1;
  else,
      h = 0;
  end;