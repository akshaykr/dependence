function [l2] = l2divergenceGine(X, Y)
% Estimates the L2 divergence between X & Y.
% TODO: kde creates a nxn matrix so don't use this for a large number of pts.

  % params
  USE_SILVERMAN_H = true;

  n = size(X, 1);
  m = size(Y, 1);
  d = size(X, 2);

  % shuffle the data
  X = X(randperm(n), :);
  Y = Y(randperm(m), :);

  n1 = round(n/2);
  n2 = n - n1;
  m1 = round(m/2);
  m2 = m - m1;
  X1 = X(1:n1, :);
  X2 = X(n1+1:end, :);
  Y1 = Y(1:m1, :);
  Y2 = Y(m1+1:end, :);

  if USE_SILVERMAN_H
    stdX = norm(std(X));
    stdY = norm(std(Y));
    hXl2 = 1.06 * stdX / n^(-1/(4+d));
    hYl2 = 1.06 * stdY / m^(-1/(4+d));
  else
    % Use Cross validation to estimate the optimal bandwidth
    [~, phat, hX] = kde(X1);
    [~, qhat, hY] = kde(Y1);
    rescale_bws = true;
    %   rescale_bws = false;
    if rescale_bws
      % Now rescale them to obtain the rates for l2 estimation
      beta = 2; % smoothness of the function 
      hXl2 = hX * n ^ (-d/( (4*beta + d) * (2*beta + d) ) );
      hYl2 = hY * m ^ (-d/( (4*beta + d) * (2*beta + d) ) );
    else
      hXl2 = hX;
      hYl2 = hY;
    end
  end

  % Now compute T1
  T1 = 0;
  for i = 1:(n1-1)
    T1 = T1 + sum( GaussKernel(hXl2, X1(i,:), X1(i+1:end, :)) );
  end
  T1 = 2 * T1/(n1 * (n1 -1)); % the h^-d term is already in GaussKernel
  % T2
  T2 = 0;
  for i = 1:(m1-1)
    T2 = T2 + sum( GaussKernel(hYl2, Y1(i,:), Y1(i+1:end, :)) );
  end
  T2 = 2 * T2/(m1 * (m1-1) );
  % T3
  G3 = GaussKernel(hXl2, X2, Y2); T3 = sum(sum(G3))/(n2*m2);
  G4 = GaussKernel(hYl2, Y2, X2); T4 = sum(sum(G4))/(n2*m2);

  l2 = sqrt( T1 + T2 - T3 - T4);

end
