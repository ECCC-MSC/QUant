function [meanx, meany] = movingmean(x,y,nsmooth)
% smooths along columns of y
for ss = 1:length(x) - nsmooth
    meanx(ss) = mean(x(ss:ss + nsmooth-1));
    meany(ss,:) = mean(y(ss:ss + nsmooth-1,:));
end

