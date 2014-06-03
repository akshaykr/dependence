function [l2] = kernel_l2(X,Y,varargin)
% Returns the l2 squared distance estimator for the two samples X,Y
% Trains a KDE on the first sample to do bandwidth selection and
% then rescales the bandwidth.
% Then uses the bi-variate kernel evaluations as the estimator.

   n1 = size(X,2);
   n2 = size(Y,2);
   d = size(X,1);

   p = inputParser;
   p.addParamValue('beta', d, @isscalar);
   p.addParamValue('debug',0, @isscalar);
   p.addParamValue('fast', 0, @isscalar);
   p.parse(varargin{:});
   Prms = p.Results;

   X1 = X(:,1:n1/2);
   X2 = X(:,(n1/2+1):end);
   Y1 = Y(:,1:n2/2);
   Y2 = Y(:,(n2/2+1):end);

   n11 = size(X1, 2);
   n12 = size(X2, 2);
   n21 = size(Y1, 2);
   n22 = size(Y2, 2);

   if Prms.fast == 0,
       %% First cross-validate the density estimate to find "optimal
       %% bandwidth for density estimation. 
       %% [est_probs, f, h_old] = kde(X1');
       

       %% We're going to throw everything away but we'll rescale the
       %% bandwidth for our problem. 
       %% This is still a hack since we don't actually know beta. 
       %% But we're trying to estimate the constant. 
       %% h = h_old*n1^(-2/(4*Prms.beta+d) + 1/(2*Prms.beta+d));
       %% h = 0.1/(n1)^(1-0.1);

       h = 0.5*n1^(-2/(4*Prms.beta+d));
       if Prms.debug == 1,
           fprintf('h_old=%0.2f h_new=%0.2f\n', h_old, h);
       end;
       
       
       T1 = GaussKernel(h, X1');
       T2 = GaussKernel(h, Y1');
       T3 = GaussKernel(h, X2', Y2');
       
       T1 = 1/(n11*(n11-1)) * sum(sum(T1 - diag(diag(T1))));
       %% fprintf('T1=%0.3f\n', T1);
       T2 = 1/(n21*(n21-1)) * sum(sum(T2 - diag(diag(T2))));
       T3 = 2/(n12*n22) * sum(sum(T3));
       if Prms.debug == 1,
           fprintf('T1=%0.3f T2=%0.3f T3=%0.3f\n', T1, T2, T3);
       end;
       l2 = T1 + T2 - T3;
   end;
   if Prms.fast == 1,
       T1 = mean(ksdensity(X1, X2));
       T2 = mean(ksdensity(Y1, Y2));
       T3 = mean(ksdensity(X2, Y2));
       l2 = T1 + T2 - 2*T3;
   end;