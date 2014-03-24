function a=join(line,s)
%  a=join(line,separator)
% Returns a string obtained by joining the elements of a cell array.
%
% See also SPLIT

% Kjartan Halvorsen
% 2001-02-22

if (nargin==1)
   s=sprintf('\t'); % tab spaced by default
end

if (isa(line,'char'))
   a=line;
   return
end

a='';

try
   a=line{1};
   for i=2:length(line)
      a=[a,s,line{i}];
   end
catch
   warning(['Unable to join: ', line])
end
