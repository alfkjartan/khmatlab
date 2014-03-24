function [m, y] = hasmissing(x)
%  [m, y] = hasmissing(x)
% Returns the number of three consecutive elements of x which are zero. This
% indicates missing marker coordinates.

% Kjartan Halvorsen
% 2003-12-03
  
m=0;

y=x;

for i=3:3:length(x)
  if ( x(i)==0 & x(i-1)==0  & x(i-2)==0 )
    m=m+1;
    y(i-2:i) = NaN;
  end
end
