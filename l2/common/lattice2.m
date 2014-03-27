function [lattice] = lattice2(interval, d),

  lattice = interval(:);

  for i = 2:d,
      T = repmat(lattice, numel(interval), 1);
      z = reshape(repmat(interval, size(lattice,1), 1), ...
                  numel(interval)*size(lattice,1), 1);
      lattice = [T z];
  end;