function [npg, gind, pind] = setc3dparam(pg, group, param)
%  npg = setc3dparam(pg, group, param)
% Searches the parameter group for the desired group and parameter
%
% Input
%    pg        ->  Parameter group struct as returned from readC3D
%    group     ->  Name of the group
%    param     ->  The parameter struct, with fields:
%                  name
%                  type
%                  data
%                  description
%                  dim
% Output
%    npg       <-  New parameter group struct

% Kjartan Halvorsen
% 2005-03-01

[slask, gind] = intersect([pg.name], group);

if isempty(gind)
  warning(['Group ', group, ' not found.']);
  p=[];
  pind = [];
  npg = pg;
  return;
end

[slask, pind] = intersect([pg(gind).Parameter.name], param.name);

if isempty(pind)
  warning(['Parameter ', param.name{1}, ' not found in group ', group, ...
	  '. Adding new parameter']);
  pind = length(pg(gind).Parameter) + 1;
end

npg = pg;

try
  npg(gind).Parameter(pind).name = param.name;
  npg(gind).Parameter(pind).datatype = param.datatype;
  npg(gind).Parameter(pind).dim = param.dim;
  npg(gind).Parameter(pind).data = param.data;
  npg(gind).Parameter(pind).description = param.description;

  
catch
  keyboard
end




