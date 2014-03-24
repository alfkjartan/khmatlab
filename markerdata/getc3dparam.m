function [p, gind, pind] = getc3dparam(pg, group, param)
%  p = getc3dparam(pg, group, param)
% Searches the parameter group for the desired group and parameter
%
% Input
%    pg        ->  Parameter group struct as returned from readC3D
%    group     ->  Name of the group
%    param     ->  Name of the parameter
% Output
%    p         <-  Empty if not found, or the parameter struct, with fields:
%                  name
%                  type
%                  data
%                  description
%                  dim

% Kjartan Halvorsen
% 2005-03-01

[slask, gind] = intersect([pg.name], group);

if isempty(gind)
  warning(['Group ', group, ' not found.']);
  p=[];
  pind = [];
  return;
end

[slask, pind] = intersect([pg(gind).Parameter.name], param);

if isempty(pind)
  warning(['Parameter ', param, ' not found in group ', group]);
  p=[];
  return;
end

p = pg(gind).Parameter(pind);




