function [mnames, np0, np0vec] = prepare_mdata(p0)
%  [mnames, np0, np0vec] = prepare_mdata(p0)
% Extracts the marker names in mdata which correspond to markers in
% p0.
%
% Input
%    p0      ->  Nested cell array with intial marker positions.
%                {{'markername1' markerdata},{...}}
% Output
%    mnames  <-  Cell array with marker names.
%    np0     <-  Nested cell array with initial marker
%                positions. The names are removed.
%                {markerdata, {...}}
%    np0vec  <-  Vectorized version of np0

% Kjartan Halvorsen
% 2003-08-04

np0 = p0;
if isempty(p0)
  np0{1} = [];
  mnames = {};
  np0vec = [];
elseif isempty(p0{1})
  mnames = {}; 
  np0vec = [];
else
  p0c=values(p0{1});
  np0{1}=cat(2,p0c{:});
  mnames = keys(p0{1});
  np0vec = cat(1, p0c{:});
end

if (length(p0) > 1)
  % recursive calls
  for br = 2:length(p0)
    [mn_br, np0_br, np0vec_br] = prepare_mdata(p0{br});
    np0{br} = np0_br;
    mnames = cat(1,mnames, mn_br);
    np0vec = cat(1, np0vec, np0vec_br);
  end
end

  