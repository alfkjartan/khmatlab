function R=rot_matrix(w,phi)
% function R=rot_matrix(w,phi)
% Returns the rotation matrix defined by the rotation of phi around 
% the (unit) vector w. If phi is not given, then it is assumed that 
% the magnitude of w is sin(phi);

% Kjartan Halvorsen
% 1999-04-22

nw=norm(w);

if (nw<1e-9)
  R=eye(3);
  return
end


if (nargin==1)
  phi = asin(norm(w));
end

w=w/nw;
W=hat(w);
R=expm_rodrigues(W,phi);
