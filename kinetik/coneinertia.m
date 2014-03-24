function [m,r0,I0,I]=coneinertia(varargin)
% [m,r0,I0,I]=inertia(height, radius1, radius2)
% Returns the moment of inertia of a truncated cone with respect to a set
% of axes at the center of the joint.
%

% Based on @twist/inertia.m
  
% Kjartan Halvorsen
% 2004-03-10
  
% The height
h = varargin{1};

% Radius of the base
r1 = varargin{2};
% Radius of the top
r2 = varargin{3};

% The cone parameter
c = r2 / r1;

% Get the mass and center of mass
[m,r0]=mass(r1,c,h);

% First the center of mass with respect to the local coordinate
% system with center at the base of the cone.
Ix=pi*r1^4*h/20*(1+c+c^2+c^3+c^4) ...
   + pi*r1^2*h^3/30*(1 + 3*c + 6*c^2);
Iz=pi*r1^4*h/10*(1+c+c^2+c^3+c^4);

I = diag([Ix; Ix; Iz]);

% Then translate the origin to the center of mass
zm=r0;
I0x=Ix-m*zm^2;
I0y=I0x;
I0z=Iz;

I0 = diag([I0x; I0y; I0z]);


function [m,rm]=mass(r,c,height)
% [m,rm]=mass(d)
% Returns the mass of the cone assuming unit density, and the
% position of the center of mass with respect to the center of the
% base of the cone. 
  
% Kjartan Halvorsen
% 2004-03-10

  % The mass
  m=pi*r^2*height*(c^2+c+1)/3;
  
  % The center of mass
  rm=height*(3*c^2+2*c+1)/(4*(c^2+c+1));
   
