function k=keys(hash)
%  k=keys(h)
% This function returns the keys of the hash 
% The hash is either a matlab 'struct' object, or 
% a 2dimensional cell array.
%
% See also values.

% Kjartan Halvorsen
% 2000-11-08

if (isa(hash,'struct'))
   k=fieldnames(hash);
elseif (isa(hash,'cell'))
   k=hash(:,1);
end
