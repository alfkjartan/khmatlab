function [tt,vv]=vect2tree(v,t)
% function [tt,vv]=vect2tree(v,t)
% Returns a tree structure equal to the structure of tt, filled
% with the contents of the vector v. Used in 'track.m'.

% Kjartan Halvorsen
% 1999-06-24

% Revisions
% 2003-09-19  Changed to work for more than one vector per level in
%             the tree

tt=t;

% Assigning the vector of the root
if (iscell(t))
  vt=t{1};
elseif (isnumeric(t))
  vt=t;
else
  error('root of tree structure must be a cell array or a numeric vector')
end

l=length(vt(:));

if (iscell(t))
  tt{1}=reshape(v(1:l), size(vt));
elseif (isnumeric(t))
  tt=reshape(v(1:l), size(vt));
end

% Stripping the vector
v(1:l)=[];
vv=v;

% recursive call
if (iscell(t))
  lt=length(t);
  if (lt>1)
    for br=2:lt
      [tt{br},vv]=vect2tree(vv,t{br});
    end
  end
end

  
