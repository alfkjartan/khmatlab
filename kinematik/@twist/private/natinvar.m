function nv=natinvar(R);
% function nv=natinvar(R);
% Computes the natural invariants of the rotationmatrix R. The
% natural invariants are a unit vector in the direction of the axis 
% and a rotation angle.

thr=1e-12;

tr=(trace(R)-1)/2;
if (tr>1)
  tr=1;
  warning('Trace of R >1')
end;

theta=acos(tr);

if(abs(theta)<thr)
  e=zeros(3,1);
else
  
  e=vect(R)/sin(theta);
end
nv=[e;theta];