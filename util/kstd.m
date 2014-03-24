function m = kstd(x, skip)
% m = kstd(x, skip)
% Computes the standard deviation of x, just as the built in function cov. The
% difference is that any rows containing only the value given in
% skip will be ignored.

% Kjartan Halvorsen
% 2003-07-31
  
if (nargin == 1 | isNan(skip))
  xx = x(~isnan(x));

else
  xx = x(find(~sum(x == skip, 2)),:)
end

m=std(xx);

