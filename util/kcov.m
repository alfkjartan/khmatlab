function m = kcov(x, y, skip)
% m = kmean(x, y, skip)
% Computes the variance of x, just as the built in function cov. The
% difference is that any rows containing only the value given in
% skip will be ignored.

% Kjartan Halvorsen
% 2003-07-31
  
if (nargin == 2)
  m = cov(x,y);
else
  xx = x(find(~sum(x == skip, 2)),:);
  m=cov(xx,y);
end
