function mmd=extractmeanmarkers(varargin)
% mmd=extractmeanmarkers(md,mnames,submnames) or
% mmd=extractmeanmarkers(mdata, subnames)
% Returns a (3 x n) matrix containing the mean marker position of
% the markers given in the cell array of marker names submnames.
% This should be a subset of the cell array mnames.
% Calls extractmarkers
  
% Kjartan Halvorsen
% 2003-11-04
%

mmd = extractmarkers(varargin{:});

mmd = mean(mmd, 1);

mmd = reshape(mmd, 3, length(mmd)/3);

