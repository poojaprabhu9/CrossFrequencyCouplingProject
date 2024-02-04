function [meandiff,Prob,meandiffs] = permTest_mean(x,y,nperm)
%permutation test to check the mean difference between x and y, and
%corresponding p- values
  x = x(:);
  y = y(:);
  lenx = length(x);
  meandiff = mean(x) - mean(y);
  z = double([x; y]);
  meandiffs = zeros(nperm,1);
  parfor j = 1:nperm
    rn = randperm(size(z,1));
    xx = z(rn(1:lenx));
    yy = z(rn(lenx+1:end));
    meandiffs(j) = mean(xx) - mean(yy);
  end
  Nbeyond = length(find(abs(meandiffs) >= abs(meandiff)));
  % Probability randomized |r| >= r_obtained
  Prob = Nbeyond/nperm;
%   disp('Done!')
%   if Prob == 0
%     disp(['p is less than' sprintf(' %.2g',1/nperm)]);
%   end
end