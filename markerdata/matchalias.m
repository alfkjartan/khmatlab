function [ma,flag]=matchalias(ali,mset1,mset2)
% [ma,flag]=matchalias(ali,mset1,mset2)
% Matches the set of names in mset1 and mset2 according to the 
% aliases in ali. Returns a hash (cell array) ma with the matching 
% names.
%
% Input
%    ali   ->  a (n x 2) cell array, where the second column is a cell
%	       array of values (aliases).
%    mset1 ->  cell array of strings.
%    mset2 ->  cell array of strings.
% Output
%    ma    <-  a hash (n x 2) cell array.

% Kjartan Halvorsen
% 2002-03-05
%
% Revisions
% 2003-01-13   the aliases in each row of ali can have different
%              lengths.
  
flag=0;

rows=size(ali,1);

ma=cell(length(mset1),2);
for k=1:length(mset1)
  for l=1:rows
    [valsfound,ind,indv]=intersect(ali{l,2},mset1{k});

    if ~isempty(valsfound)

      for j=1:length(mset2)
	[matchfound,indarow]=intersect(ali{l,2},mset2{j});
	if ~isempty(matchfound)
	  ma{k,1}=mset1{k};
	  ma{k,2}=mset2{j};
	  break
	end
      end
    
    end
  end
  if isempty(ma{k,1})
    msg=['Marker  ', mset1{k},...
	 '  not found in alias table.'];
    disp(msg);
    flag=1;
  end
end



