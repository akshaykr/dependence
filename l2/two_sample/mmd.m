function [h,p,that] = mmd(X,Y),
    addpath('../common/');
    addpath('./mmd/');
    %% unbiased MMD with gaussian kernel with fixed bandwidth (fixed at 1).
    %% Using Arthur's code
    params = [];
    params.sig = 1;
    params.numEigs = -1;
    params.numNullSamp = size(X,2);

    [that, thresh] = mmdTestSpec(X',Y',0.05, params);
    if that > thresh,
        h = 1;
    else,
        h = 0;
    end;

    p = 0;