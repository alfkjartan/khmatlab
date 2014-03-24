function [v,th]=rotvect(R)
% function v=rotvect(R)
% The rotation vector associated with a rotation matrix.

% Kjartan Halvorsen
% 2007-01-11

tolSin = 1 + 1e-10;

RI = R-eye(3);
[V,D] = eig(RI);

for ii=1:3
  if isreal(D(ii,ii))
    ind = ii;
    break
  end
end

v = V(:,ind);
restind = setdiff(1:3, ind);
v2 = real(sum(V(:,restind),2));
v2 = v2 / norm(v2);

sinth = (cross(v2,R*v2))'*v;
if abs(sinth) > tolSin 
  th = pi/2;
  v = sign(sinth)*v;
  return
else
  th = asin(sinth);
end

if th<0
  th = -th;
  v = -v;
end

