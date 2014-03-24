function [hh,aliasfound]=reducealias2(h,vals)
%  hh=reducealias2(h,vals)
% Removes values associated with keys that are not contained in vals.
% Modified version of reducehash. Modified again 2002-03-24 and 2002-12-16
%
% Input
%    h     ->  a (n x 2) cell array, where the second column is a cell
%	       array of values (aliases).
%    vals  ->  cell array of strings, or possibly cell array of cell
%	       arrays. 
% Output
%    hh    <-  a hash (n x 2) cell array.

% Kjartan Halvorsen
% 2002-03-24

rows=size(h,1);

hh={};

aliasfound=zeros(rows,1); 
    % Will at the end contain a 1 for each row where an alias is
    % found.

for k=1:rows
   [valsfound,ind,indv]=intersect(h{k,2},vals);
   ind=sort(ind);
   indv=sort(indv);

   if isempty(valsfound)
     %      msg=['Alias for marker ', sprintf('%s  ',h{k,1}), ...
     %           'not found in data. ', sprintf('\n'),...
     %           ' Problems may occur later.',sprintf('\n')];
     %      warning(msg);
     hh{k,1}=h{k,1};
     hh{k,2}=h{k,1};
      
   else
      hh{k,1}=h{k,1};
      hh{k,2}=h{k,2}{ind(1)};
      aliasfound(k)=1;
   end
end



