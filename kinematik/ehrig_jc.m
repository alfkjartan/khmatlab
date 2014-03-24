function [c1, c2] = ehrig_jc(p1, p2)
%% [c1, c2] = ehrig_jc(p1[, p2])
%% Implementation of the SCoRe method by Ehrig et al 2006.
%%
%% Input
%%   p1    ->  data set Nfr x Nmarkers1*3 for segment 1
%%   p2    ->  data set Nfr x Nmarkers2*3 for segment 2. Optional if not given, the center of
%%             rotation is fixed in the lab coordinate system.
%%
%% Output
%%   ci     <-  the center of rotation in local coordinate system i. That is, as local
%%              coordinate system which coincides with the global system for the first 
%%              frame of data.

if nargin == 0
   do_unit_test();
else

  [Nfrs, Nmrks3] = size(p1);

  if nargin == 1
    T2 = zeros(Nfrs, 16);
    T2(:, [1 6 11 16]) = 1;
  else
    T2 = getMotion(p2);
  end

  T1 = getMotion(p1);

  A = zeros((Nfrs-1)*3, 6);
  b = zeros((Nfrs-1)*3, 1);

  for i=2:Nfrs
    TT1 = reshape(T1(i,:), 4, 4);
    TT2 = reshape(T2(i,:), 4, 4);

    A( (i-1)*3+1:i*3, 1:3 ) = TT1(1:3, 1:3);
    A( (i-1)*3+1:i*3, 4:6 ) = -TT2(1:3, 1:3);

    b( (i-1)*3+1:i*3) = TT2(1:3, 4) - TT1(1:3, 4);
  end

  cc = A \ b;

  c1 = cc(1:3);
  c2 = cc(4:6);
end

function do_unit_test()

%% Create som test data. Points on rigid body moving about the point (1,1,1);
c_true = [1;1;1];

Nfrs = 20;
ampl = 30*pi/180;

phi = linspace(0,ampl, Nfrs);
th = linspace(0, ampl/2, Nfrs);


%% Twist 1
w1 = [1;0;0]; % The direction
v1 = -cross(w1, c_true);
tw1 = hat([v1;w1]);

%% Twist 2
w2 = [0;1;0]; % The direction
v2 = -cross(w2, c_true);
tw2 = hat([v2;w2]);

% Points
p4 = [2 1 0 0
      0.5 2 2 0
      0 0 2 2
      1 1 1 1];

%% Move points
Nmarks = size(p4,2);
md = zeros(Nfrs, Nmarks*3);

for i=1:Nfrs
    g1 = expr(tw1*phi(i));
    gg = g1*expr(tw2*th(i));

    pp = gg*p4;

    md(i,:) = reshape(pp(1:3, :), 1, Nmarks*3);
end

c = ehrig_jc(md);

disp("Unit test 1");
if max(c-c_true) > 1e-12
  disp("Failed")
  disp("Expected"); c_true
  disp("Found"); c
else
    disp("OK")
end

     
     
 
