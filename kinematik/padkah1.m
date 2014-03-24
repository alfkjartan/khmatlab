function [th,flag]=padkah1(tw,p,q)
% function [th,flag]=padkah1(tw,p,q)
% Computes the solution (if any) to the Paden-Kahan subproblem
% 1. See Murray, Li, Sastry "Robotic Manipulation", 99-100
% th is the angle of rotation around the given twist that rotates
% the point p into q.
%
% Input
%    tw      ->   twist
%    p       ->   a point
%    q       ->   the same point after the rotation.
%
% Output
%    th      <-   the angle of rotation around the twist axis.
%    flag    <-   1 if no solution
%

% Kjartan Halvorsen
% 1999-09-21
% Modified 2003-09-17, to be used with simple twist matrices

tol1=1e-8;

v=tw(1:3,4);
w=vect(tw(1:3,1:3));

% a point on the axis
r=cross(w,v);

u=p-r;
v=q-r;

% projector onto the plane normal to w:
pii=eye(3)-w*w';

% project u,v onto the plane
u_p=pii*u;
v_p=pii*v;

% check for existence of solution
tol2=1e-3;
if (abs(norm(u_p)-norm(v_p))>tol2)
  flag=1;
elseif (abs(w'*u-w'*v)>tol2)
  flag=1;
else
  flag=0;
end

uxv=cross(u_p,v_p);
w_uxv=w'*uxv;
uv=u_p'*v_p;

if (abs(uv)<tol1)
  th=asin(w_uxv);
else
  th=atan(w_uxv/uv);

  % Resolve ambiguity of atan in [-pi,pi]
  if (sign(w_uxv)==-1 & sign(uv)==-1)  %-pi<th<-pi/2
    th=th-pi;
  elseif (sign(w_uxv)==1 & sign(uv)==-1) % pi/2<th<-pi
    th=th+pi;
  end
  
  if (abs(th)<tol1)
    if (uv<0)
      th=pi;
    end
  end
end

  
