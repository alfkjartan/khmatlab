function ismd = ismarkerdata(md)
%  ismd = ismarkerdata(md)
% Returns true if md is a marker data struct: {attr, markerdata}.
  
% Kjartan Halvorsen
% 2004-01-27
  
ismd = (iscell(md) & iscell(md{1}) & ...
	(size(md{1},2)==2) & isnumeric(md{2}));

