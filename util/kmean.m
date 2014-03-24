function m = kmean(x, skip)
% m = kmean(x, skip)
% Computes the mean of x, just as the built in function mean. The
% difference is that any rows containing only the value given in
% skip will be ignored.

% Kjartan Halvorsen
% 2003-07-31
  
if (nargin == 1 | isnan(skip))
  xx = x(~isnan(x));

else
  xx = x(find(~sum(x == skip, 2)),:)
end

m=mean(xx);



