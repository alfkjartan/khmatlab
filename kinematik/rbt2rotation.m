function m=rbt2rotation(rbt)
% function m=rbt3vrmlform(rbt)
% Extracts the rotation vector and angle, returned in 1x4 row vector,
% from rbt.

% Kjartan Halvorsen
% 2000-05-03

if (size(rbt,2)==4)
  rot=rbt(1:3,1:3);
else
  rot=rbt;
end

rasym=0.5*(rot-rot');
rotv=vect(rasym);
rn=norm(rotv);
if (rn>1e-12) rotv=rotv/rn; end
  
rotv=rotv(:);

tr=trace(rot);
if (tr>3) tr=3; end
theta=acos(0.5*(tr-1));

m=[rotv' theta];
