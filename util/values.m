function c=values(h)
%  c=values(h)
% This function returns the data in the hash as a cell array.
% The hash is either a matlab 'struct' object, or 
% a 2dimensional cell array.
%
% See also keys

% Kjartan Halvorsen
% 2000-11-24

if (isa(h,'struct'))
   c=struct2cell(h);
elseif (isa(h,'cell'))
   c=h(:,2);
end
