function [r, pred, drdx]=residual_rigidbody(y, x, q, dt, d0, r0)
%%  [r, pred, drdx]=residual_rigidbody(y, x, q, dt, d0, r0)
%%
%% Measurement function for a rigid body:
%%   y1 = d0 + d + q*exp(dt*w)*r1*exp(-dt*w)*qconj
%%
%% Input
%%    y      ->   the current measurement y = [y1;y2; ...;ym]
%%    x      ->   the current state. x = [w; alpha; d; v; a]
%%    q      ->   the previous orientation. The current orientation is
%%                  qn = q*exp(dt*w).
%%    dt     ->   the sampling time
%%    d0     ->   the reference position of the centroid (3x1) 
%%    r0     ->   the reference position of the marker vectors (from
%%                centroid to marker) (4 x m) 
% Output
%    r       <-   the residual y - pred

% Kjartan Halvorsen
% 2012-05-29
%
%

if (nargin == 0)
  do_unit_test();
else
  
  w = x(1:3);
  d = x(7:9);

  qc = qconj(q);
  pred = y;

  nmarkers = size(r0,2);

  drdw = zeros(nmarkers*3, 3);

  [ew, dewdw] = qexp(dt*w);
  qe = qmult(q, ew);
  qec = qconj(qe);
  
  dqdw = zeros(4,3);
  for i=1:3
    dqdw(:,i) = qmult(q, dewdw(:,i));
  end


  for m=1:nmarkers
    [rm, drmdw] = qrotexp(r0(:,m), q, w, qc, qe, qec, dqdw);
    pred((m-1)*3+1:m*3) = d0 + d + rm(1:3);
    drdw((m-1)*3+1:m*3,:) = -dt*drmdw(1:3,:);
  end

  r = y-pred;

  drdx = zeros(3*nmarkers, length(x));
  drdx(:,1:3) = drdw;
  drdx(:,7:9) = -repmat(eye(3), nmarkers, 1);
 
end

function do_unit_test()
  disp("Unit test for function residual_rigidbody")

  tol = 1e-6;
  dx = 1e-8;
  dt = 0.01;
  dt = 1;

  q = (quaternion(randn(3,1), rand(1)))';

  %x = 10*pi/180*randn(9,1);
  x  = zeros(9,1);
  x(4:end) = randn(6,1);
  x(7:9) = 200*pi/180*randn(3,1);
  I9 = eye(9);

  p0 = randn(3,4);
  d0 = mean(p0,2);
  r0 = cat(1, p0 - repmat(d0, 1, 4), zeros(1,4));

  qe = qpropagate(q, x(1:3), dt);
  y = zeros(3,4);
  for m=1:4
    y(:,m) = d0 + x(4:6) + qtransv(r0(1:3,m), qe);
  end
  y = y(:);

  [r, pred, drdx] = residual_rigidbody(y, x, q, dt, d0, r0);
  
  dr_dx = zeros(12,9);
  for i = 1:9
    dr_dx(:,i) = ( residual_rigidbody(y, x+dx*I9(:,i), q, dt, d0, r0) - r ) / dx;
  end
  
  if ( norm(drdx - dr_dx) > tol )
    disp('Test 1. Failed')
    cat(2, drdx, dr_dx)
    keyboard
  else
    disp('Test 1. OK')
  end

