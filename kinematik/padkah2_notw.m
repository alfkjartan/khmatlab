function [th,flag,thalt]=padkah2_notw(r,w1,w2,p,q)
% function[th,flag,thalt]=padkah2_notw(r,w1,w2,p,q)
% Computes the solution (if any) to the Paden-Kahan subproblem
% 2. See Murray, Li, Sastry "Robotic Manipulation", 100-102.
% th=[th1;th2] are the two rotations around the two twists
% intersecting at r and with directions w1 and w2.
% Based on @twist/padkah2.
% Input
%    r       ->   intersection of twist axes.  
%    w1,w2   ->   unit vectors in direction of twist axes
%    p       ->   a point
%    q       ->   the same point after the rotation.
%
% Output
%    th      <-   the two angle of rotation.
%    flag    <-   1 if no solution, 0 if ok.
%    thalt   <-   the alternative solution.
%

% Kjartan Halvorsen
% 2004-03-19
  
% The outline of the solution is to find a point c which is such
% that p is first rotated an angle th2 around twist2 to c, and then 
% c is rotated an angle th1 around twist1 to the point q.
% Once c is found, use Paden-Kahan subproblem 1 ('twist/padkah1')
% to find the solution.

u=p-r;
v=q-r;

alpha = ((w1'*w2)*(w2'*u) - w1'*v) / ((w1'*w2)^2 - 1);
beta = ((w1'*w2)*(w1'*v) - w2'*u) / ((w1'*w2)^2 - 1);
gamma2 = ((norm(u))^2 - alpha^2 - beta^2 -2*alpha*beta*w1'*w2) / ...
    (norm(cross(w1,w2)))^2;
gamma11=sqrt(gamma2);
if (~isreal(gamma11))
  warning('No solution to the paden-kahan subproblem type 2')
  th=[0;0]
  flag=1;
  thalt = [NaN;NaN];
else
  gamma12=-gamma11;

  z1 = alpha*w1 + beta*w2 + gamma11*cross(w1,w2);
  z2= alpha*w1 + beta*w2 + gamma12*cross(w1,w2);

  c1 = z1 + r;
  c2 = z2 + r;

  % Now use padkah1
  [th21,flag21]=padkah1_notw(r,w2,p,c1);
  [th11,flag11]=padkah1_notw(r,w1,q,c1);

  [th22,flag22]=padkah1_notw(r,w2,p,c2);
  [th12,flag12]=padkah1_notw(r,w1,q,c2);

  % By definition choose the solution with smallest combined
  % magnitude of rotations
  if (abs(th21) + abs(th11) < abs(th22) + abs(th12))
    th=[-th11;th21];
    thalt=[-th12;th22];
    flag= flag11 | flag21;
  else
    th=[-th12;th22];
    thalt=[-th11;th21];
    flag= flag12 | flag22;
  end
end

