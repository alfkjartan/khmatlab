function y = rotate_spherical_coords(x,R) 
%  y = rotate_spherical_coords(x,R) 
% Rotates the unit vector parameterized by spherical coordinates.
%
% Input
%    x    -> spherical coordinates x=[th;ph], so that the unit
%            vector is w = [sin(ph)*cos(th); sin(ph)*sin(th);
%                           cos(th)]
%    R    -> Rotation matrix
% Output
%    y    -> spherical coordinates y=[th;ph]

% Kjartan Halvorsen
% 2007-06-18

if (nargin > 0)
  w = [sin(x(2))*cos(x(1)); sin(x(1))*sin(x(2)); cos(x(2))];

  w2 = R*w;

  y = [atan2(w2(2), w2(1)); acos(w2(3))];

else 
  % unit test
  
  x = [0;pi/2];
  R = expm_rodrigues(hat([0;0;1]), pi/2);
  
  y_expected = [pi/2;pi/2]
  y_found = rotate_spherical_coords(x,R)
  
end
