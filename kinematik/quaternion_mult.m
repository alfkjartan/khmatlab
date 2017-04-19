function q3 = quaternion_mult(q1,q2)
%  q3 = quaternion_mult(q1,q2)
% Quaternion (Grassmann) multiplication
% 
% Input
%    qi   ->  quaternion. Either (4x1) real vector or (2x2) imaginary matrix
% Output
%    q3   <-  quaternion. Either (4x1) real vector or (2x2) imaginary matrix

% Kjartan Halvorsen
% 2007-06-18

if (nargin > 0)
  if (length(q1) == 4)
    if 0
      q1 = [q1(1)+i*q1(2) q1(3)+i*q1(4)
	       -q1(3)+i*q1(4) q1(1)-i*q1(2)];
      q2 = [q2(1)+i*q2(2) q2(3)+i*q2(4)
	      -q2(3)+i*q2(4) q2(1)-i*q2(2)];
      q3 = q1*q2;
    
      q3r = real(q3);
      q3i = imag(q3);
      q3 = [q3r(1,1); q3i(1,1); q3r(1,2); q3i(1,2)];
    end
 
    
    r1 = q1(4);  % scalar part
    v1 = q1(1:3); % vector part
    r2 = q2(4);
    v2 = q2(1:3);
 
    v1 = v1(:);
    v2 = v2(:);
    
    r3 = r1*r2 - dot(v1,v2);
    v3 = r1*v2 + r2*v1 + cross(v1,v2);
    
    q3 = cat(1, v3(:), r3);
  else
    q3 = q1*q2;
  end

else
  % Unit test

  twobytwo = 0;
  
  w1 = [1;0;0];
  w2 = [0;1;0];
  R1 = expm_rodrigues(hat(w1), pi/3);
  R2 = expm_rodrigues(hat(w2), pi/5);
  
  
  q1 = rotation2quaternion(R1,twobytwo);
  q2 = rotation2quaternion(R2,twobytwo);
  
  q3 = quaternion_mult(q1,q2)
  
  q3expected = rotation2quaternion(R1*R2, twobytwo)
  
  
  
  
end
