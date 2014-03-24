function twn=setTheta(tw,th)
% function twn=getTheta(tw,th)
% Sets the rotation angle of the twist (screw)

% Kjartan Halvorsen
% 1999-06-01

tw.theta=th;
twn=twist(tw);