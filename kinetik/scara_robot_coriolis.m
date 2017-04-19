function C = scara_robot_coriolis(sm, th, thdot)
%  C = scara_robot_coriolis(sm, th, thdot)
% Hardcoded from p178, Lee, Murrray, Sastry.

M1 = sm.inertia{1};
M2 = sm.inertia{2}{1};
M3 = sm.inertia{2}{2}{1};
M4 = sm.inertia{2}{2}{2}{1};

m1 = sm.mass{1};
m2 = sm.mass{2};
m3 = sm.mass{3};
m4 = sm.mass{4};

l1 = sm.l1;
l2 = sm.l2;

%% From p 178
gamma = l1*l2*m3 + l1*l2*m4 + l1*l2/2*m2;
s2 = sin(th(2));

C = zeros(4,4);
C(1,1) = -gamma*s2*thdot(2);
C(1,2) = -gamma*s2*(thdot(1)+thdot(2));
C(2,1) = gamma*s2*thdot(1)

