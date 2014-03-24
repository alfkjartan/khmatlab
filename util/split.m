function a=split(line,s)
%  a=split(line,separator)
% Returns a cell array of strings obtained by splitting line using the
% specified separator.
%
% See also JOIN

% Kjartan Halvorsen
% 2001-02-08

if (nargin==1)
   s=sprintf('\t'); % tab spaced by default
end

[tok,rem]=strtok(line,s);

a{1}=tok;

k=1;
while (~isempty(rem))
   k=k+1;
   rem=deblank(rem);
   rem=deblank(fliplr(rem));
   rem=fliplr(rem);
   [tok,rem]=strtok(rem,s);
   a{k}=tok;
end

