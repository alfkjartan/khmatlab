function shaft_vibration(EI, a, b, c, m)

if (nargin == 0)
  do_unit_test();
else
  L = (2*a+b);
  my = m/L; % unit kg/m
  L0 = b;
  
  omega = sqrt(EI/my)*pi^2/L0^2 % Frequency of first mode of vibration

  dt = 2*pi/omega
end




function w = beam_deflection(x, L, c1, c2, c3, c4)
  xL = x*L;
  w = c1*sin(xL) + c2*cos(x*L) + c3*sinh(xL) + c4*cosh(xL);


function do_unit_test()
  a = 0.5;
  b = 0.6;
  c = 0.075;
  m = 1.0;

  %% EI from geometry
  r = 15e-3; % radius of the shaft
  t = 15e-4; % thickness of wall
  II = pi*r^3*t; %% Thin tube approximation to second area of moment
  E = 240e9; % in GPa = 10^9 N / m^2
  EI = E*II; % Unit N*m^2 or kg m/s^2 m^2
  
  shaft_vibration(EI, a,  b, c, m);


    