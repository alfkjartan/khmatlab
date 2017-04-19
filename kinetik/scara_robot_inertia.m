function [M, M1, M2, M3, M4] = scara_robot_inertia(sm, th2)
%  M = scara_robot_inertia(sm, th2)
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

if length(th2) > 1
    th2 = th2(2);
end

%% From p 177
alpha = M1(6,6) + (l1/2)^2*m1 + l1^2*(m2+m3+m4);
beta = M2(6,6) + M3(6,6) + M4(6,6) + (l2/2)^2*m2 + l2^2*(m3+m4);
gamma = l1*l2*m3 + l1*l2*m4 + l1*l2/2*m2;
delta = M3(6,6) + M4(6,6);
 
M = [alpha+beta+2*gamma*cos(th2) beta+gamma*cos(th2) delta 0
     beta+gamma*cos(th2) beta  delta 0
     delta delta delta 0
     0  0  0 m4];

