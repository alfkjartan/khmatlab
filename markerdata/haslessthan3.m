function [m,xx] = haslessthan3(x)
%  [m,xx] = haslessthan3(x)
% Returns 1 if the set of markers has less than three non-missing markers.

% Kjartan Halvorsen
% 2003-02-04

[mssng, xx] = hasmissing(x);
m = ( (length(x)/3 - mssng) < 3 ) ;
