function v=tree2vect(t)
% function v=tree2vect(t)
% Returns a column vector obtained by stacking the contents of the
% tree. Used in 'track.m'.

% Kjartan Halvorsen
% 1999-06-24

v=t{1};
v=v(:);

lt=length(t);
if (lt>1)
  for br=2:lt
    vb=tree2vect(t{br});
    v=cat(1,v,vb);
  end
end

  
