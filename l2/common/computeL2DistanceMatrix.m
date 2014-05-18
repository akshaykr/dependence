function [D] = computeL2DistanceMatrix(data)
% data is a cell array containing num_dists matrices. Each matrix of size
% nixnum_dims gives the points sampled from each distribution.

  % Prelims
  num_dists = numel(data);
  
  % First compute the similarity matrix
  D = zeros(num_dists);
  for i = 1:num_dists
    if (num_dists > 40) && (mod(i,10) == 0)
      % then report progress
      fprintf('Computing L2 between %dth distribution and others\n', i);
    end
    for j = (i+1):num_dists
%       l2ij = l2Divergence(data{i}, data{j});
      l2ij = l2DivergenceGine(data{i}, data{j});
      if imag(l2ij) > 0, l2ij = 0;
      end
      D(i,j) = l2ij;
      D(j,i) = l2ij;
    end
  end
end
