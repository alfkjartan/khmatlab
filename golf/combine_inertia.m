function [II, CoM, mass]  = combine_inertia(I1, CoM1, mass1, I2, ...
                                            CoM2, mass2, refPoint)
%%  [II, CoM, mass]  = combine_inertia(I1, CoM1, mass1, I2, CoM2, mass2)
% Combines the inertial parameters of two parts of the same rigid body. The CoM and
% inertia matrices must be with respect to the same coordinate system.
%
% If refPoint is not given, compute wrt combined CoM.

%% Kjartan Halvorsen
% 2013-08-23

if (nargin == 0)
   do_unit_test();
   return;
end

mass = mass1 + mass2;

if mass < 1e-12
   CoM = zeros(3,1);
else
  CoM = (CoM1*mass1 + CoM2*mass2) / mass;
end

if nargin < 7 
 v1 = CoM - CoM1;
 v2 = CoM - CoM2;
else
    v1 = refPoint - CoM1;
    v2 = refPoint - CoM2;
end

%%II = I1 + I2 + mass1 * diag( [v1(2:3)'*v1(2:3)
%			      v1([1 3])'*v1([1 3])
%			      v1(1:2)'*v1(1:2)] ) ...
%    + mass2 * diag( [v2(2:3)'*v2(2:3)
%			      v2([1 3])'*v2([1 3])
%			      v2(1:2)'*v2(1:2)] );
II = I1 + I2 + mass1 * parallel_axis(v1) ...
             + mass2 * parallel_axis(v2);


function JJ = parallel_axis(t)
%% The unscaled change in moment of inertia because of translation of reference point
% t is the vector from old point to new
%JJ = [-(t(2)^2+t(3)^2) t(1)*t(2) t(1)*t(3)
%      t(1)*t(2) -(t(1)^2 + t(3)^2) t(2)*t(3)
%      t(1)*t(3) t(2)*t(3) -(t(1)^2 + t(2)^2)];
 
JJ = t'*t*eye(3) - t*t';

function do_unit_test()

%% Create two bodies. Known geometries
ex = randn(3,1);
ex = ex / norm(ex);

ey = randn(3,1);
ey = ey - (ey'*ex)*ex;
ey = ey/norm(ey);

ez = cross(ex, ey);

m1 = 1;
m2 = 2;

com1 = [1;0;0];
com2 = [-1;0;0];

I1 = eye(3);
I1(1,3) = 1.1111;
I1(3,1) = 1.1111;

I2 = eye(3);
I2(1,2) = 1.222;
I2(2,1) = 1.222;

[II, com, mass] = combine_inertia(I1, com1, m1, I2, com2, m2);

assert(mass, m1+m2, 1e-14);
assert(com, [-1/3;0;0], 1e-14);

II


