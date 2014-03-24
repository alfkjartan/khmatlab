function [m,xx] = haslessthann(x,n)
%  [m,xx] = haslessthann(x,n)
% Returns 1 if the set of markers has less than n non-missing markers.
% Generalization of haslessthan3.

% Kjartan Halvorsen
% 2005-02-28

[mssng, xx] = hasmissing(x);
m = ( (length(x)/3 - mssng) < n ) ;
