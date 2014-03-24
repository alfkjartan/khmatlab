function q = rot2quat(R)
%%  q = rot2mat(R)
%% Creates a quaternion from a rotation matrix

%% 2012-03-27
%% Kjartan Halvorsen

if (nargin == 0)
  do_unit_test();
else
  tol = 1e-8;

  t = trace(R);
  r = sqrt(1+t);
  s = 0.5/r;
  w = 0.5*r;

  x = (R(3,2)-R(2,3))*s;
  y = (R(1,3)-R(3,1))*s;
  z = (R(2,1)-R(1,2))*s;

  q = cat(1, x,y,z, w);
end

function do_unit_test()
  disp('Unit test of rot2quat')

  ex = randn(3,1);
  ex = ex / norm(ex);
  ey = randn(3,1);
  ey = ey - (ey'*ex)*ex;
  ey = ey / norm(ey);
  ez = cross(ex,ey);

  R = cat(2, ex, ey, ez);
  
  q = rot2quat(R);
  R1 = qtransvmat(q);

  if (norm(R-R1) > 1e-12)
    disp('Test 1: Failed')
    disp('Expected'), disp(R)
    disp('Found'), disp(R1)
  else
    disp('Test 1: OK')
  end
