function [p,flag]=intersection(tw1,tw2)
% function [p,flag]=intersection(tw1,tw2)
% Computes the point of intersection of two twists. If the twists
% do not intersect, then the point closest to the two lines is
% returned and the flag set to 1.

% Kjartan Halvorsen
% 1999-09-21
% Modified 2003-09-17 to work for simplw twist matrices
  
v1=tw1(1:3,4);
w1=vect(tw1(1:3,1:3));
q1=cross(w1,v1);

v2=tw2(1:3,4);
w2=vect(tw2(1:3,1:3));
q2=cross(w2,v2);


% A line (as the axis of the twist) can be described as the
% intersection of two planes. The point of intersection of the two
% axes is consequently the intersection of four planes.

[U,S,V]=svd(w1);
z1=U(:,find(S==0));
[U,S,V]=svd(w2);
z2=U(:,find(S==0));

A=[z1';z2'];
b=[z1'*q1;z2'*q2];

p=A\b; % The point of intersection.

tol=1e-9;
if (norm(A*p-b)>tol)
  flag=1; % The twists don't intersect
else
  flag=0;
end

