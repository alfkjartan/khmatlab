function [th, flag] = invkin_twoaxes(p, q, w1, w2)
%  [th, flag] = invkin_twoaxes(p, q, w1, w2)
% 
% Solves 
%     q= R1*R2*p = exp{hat(w1)*th1} * exp{hat(w2)*th2) * p,
% 
% input
%    p     ->  (3 x m) set of markers
%    q     ->  (3 x m) set of markers
%    w1    ->  (3 x 1) unit vector in the direction of the first axis.
%    w2    ->  (3 x 1) unit vector in the direction of the second
%              axis. The axes must be perpendicular.
% output
%    th    <-  angles [th1;th2]
%    flag  <-  1 indicates solution not found
  
% Kjartan Halvorsen
% 2004-03-18

tol = 1e-3;
pn = norm(p);
qn = norm(q);

if (abs(pn-qn) > tol)
  th=[0;0];
  flag = 1;
  return
end

% th2 is angle between q and the plane normal to w1

% Project q onto plane normal to w1.

qpr = (eye(3) - w1*w1')*q;

th2 = sign(w2'*cross(qpr,q))*acos(qpr'*q/(qn*norm(qpr)));

% th1 is angle between qpr and p
th1 = sign(w1'*cross(p,qpr))*acos(qpr'*p/(pn*norm(qpr)));

th = [th1;th2];
flag = 0;
