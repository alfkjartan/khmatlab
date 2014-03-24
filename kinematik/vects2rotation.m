function R=vects2rotation(e1,e2)
% function R=vects2rotation(e1,e2)
% Returns the rotation matrix that rotates e1 into e2.

% Kjartan Halvorsen
% 1999-04-22

e1=e1(:);
e2=e2(:);

e1=e1/norm(e1);
e2=e2/norm(e2);

w=cross(e1,e2);
theta=acos(e1'*e2);

if(norm(w)<1e-8)
  if (sign(e1'*e2)==1)
    R=eye(3);
  else % The vectors are in opposite directions. Rotate pi radians
    % around a vector normal to e1.
    ee1=e1;
    II=eye(3);
    ee=II(:,find(e1==min(abs(e1))));
    ee=ee(:);
    ee=ee(1:3);
    ee2=ee - (ee'*ee1)*ee1;
    ee2=ee2/norm(ee2);
    ee3=cross(ee1,ee2);
    % Base 1:
    E1=[ee1 ee2 ee3];
    % Base 2:
    E2=[-ee1 ee2 -ee3];
    % Seek R, R*E1=E2 => R=E2*E1'
    R=E2*E1';
  end
else
  R=rot_matrix(w,theta);
end
