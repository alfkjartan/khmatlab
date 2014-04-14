function [alpha, sep] = ellipsoid_direction(W, v)
%%  [alpha, sep] = ellipsoid_direction(W, v)
%% Calculates the direction between the main axis of the mobility ellipsoid and the vector v.
%%
%% Input
%%   W      ->  Mobility matrix 3 x 3 
%%   v      ->  Vector. Need not be unit length
%% Output
%%   alpha  <-  Angle in degrees
%%   sep    <-  Measure of how elongated the ellipsoid is. It is the relative separation of 
%%              the first and second singular value of W: (sigma_2-sigma_1)/sigma_1. For a 
%%              ellipsoid collapsed to one dimension, this value is 1. For a completely round,
%%              or squeezed (round but flat) ellipsoid the value is 0.

%% Kjartan Halvorsen
%% 2014-04-14

if nargin==0
   do_unit_test();
else
  [U,S,V] = svd(W);

  sep = (S(1,1) - S(2,2)) / S(1,1);

  ev = v(:) / norm(v);

  alpha = acosd(V(:,1)'*ev);

  if alpha > 90 
    alpha = 180-alpha;
  end
end

function do_unit_test()

W = eye(3)
W(1,1) = 2;

%% Rotation about z-axis
th = 42;
R = [cosd(th) sind(th) 0
     -sind(th) cosd(th) 0
     0  0 1];

WW = R'*W*R;

[alpha,sep] = ellipsoid_direction(WW,[3;0;0])

assert(norm(alpha-th)<1e-12)
disp('Test 1 OK!')

assert(norm(sep-0.5)<1e-12)
disp('Test 2 OK!')



