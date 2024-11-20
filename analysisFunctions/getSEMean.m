% [se,bs] = getSEMean(X,N) returns the standard error of the mean of
% the values in X. 

function [se,bs] = getSEMean(X,N)

if ~exist('N','var');                   N=length(X);                    end

bs = bootstrp(N,@mean,X);
se = std(bs);
end