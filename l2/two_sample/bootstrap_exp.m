% Script for testing bootstrap
addpath('../common/');

B = 100;

theta = 1.0;
C100 = zeros(100,2);
c100 = 0;
C200 = zeros(100,2);
c200 = 0;
C400 = zeros(100,2);
c400 = 0;
%% 100 samples
for i=1:100,
    X = rand([1,500]);
    Y = rand([1,200]);
    [that,vals] = bootstrap_square(X,B);
    s = sort(real(vals-that));
    C100(i,1) = that - s(int64(0.95*B));
    C100(i,2) = that - s(int64(0.05*B));
    if theta >= C100(i,1) && theta <= C100(i,2),
        c100 = c100+1;
    end;
    fprintf('100: iter: %d that: %0.2f, interval: (%0.2f, %0.2f)\n', ...
            i, that, C100(i,1), C100(i,2));
end;
