function [r, pred, drdx]=observe_box(y, x, q, dt, d0, r0)
%%  [r, pred, drdx]=observe_box(y, x, q, dt, d0, r0)
%%
%% Measurement function for the box (a rigid body):
%%   y1 = d0 + d + q*exp(dt*w)*r1*exp(-dt*w)*qconj
%%
%% Input
%%    y      ->   the current measurement y = [y1;y2; ...;ym]
%%    x      ->   the current state. x = [w; d; v]
%%    q      ->   the previous orientation. The current orientation is
%%                  qn = q*exp(dt*w).
%%    dt     ->   the sampling time
%%    d0     ->   the reference position of the centroid (3x1) 
%%    r0     ->   the reference position of the marker vectors (from
%%                centroid to marker) (4 x m) 
% Output
%    r       <-   the residual y - pred

% Kjartan Halvorsen
% 2012-03-20
%
%

if (nargin == 0)
  do_unit_test();
else
  
  w = x(1:3);
  d = x(4:6);

  %%[ew, dedw, decdw] = qexp(w*dt);
  %%[ewc, decdw] = qexp(-w*dt);
  
  qc = qconj(q);
  %%qe = qmult(q,ew);
  %%qec = qconj(qe);

  pred = y;
  %r = y-pred;

  %%qdedw = dedw;
  %%decdwqc = decdw;
  %%for i=1:3
  %%  qdedw(:,i) = qmult(q, dedw(:,i));
  %%  decdwqc(:,i) = qmult(decdw(:,i), qc);
  %%end

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
    %%rm = qmult(qe, qmult(r0(:,m), qec));
    %%[rm, drmdq] = qrot1(r0(:,m), qe, qec);
    [rm, drmdw] = qrotexp(r0(:,m), q, w, qc, qe, qec, dqdw);
    
    pred((m-1)*3+1:m*3) = d0 + d + rm(1:3);
    
    %%drmdw = zeros(4,3);
    %%for i=1:3
    %% drmdw(:,i) = -dt*(qmult(qdedw(:,i), qmult(r0(:,m), qc)) - \
    %%  qmult(qmult(q, r0(:,m)), decdwqc(:,i)));
    %%end
    %%ddw = dt*drmdq*dedw;
    %%drdw((m-1)*3+1:m*3,:) = ddw(1:3,:);

    drdw((m-1)*3+1:m*3,:) = -dt*drmdw(1:3,:);
  end

  r = y-pred;

  drdx = cat(2,  drdw, \
	     -repmat(eye(3), nmarkers, 1), \
	     zeros(nmarkers*3, 3));
 
end

function do_unit_test()
  disp("Unit test for function observe_box")

  tol = 1e-6;
  dx = 1e-8;
  dt = 0.01;
  dt = 1;

  q = (quaternion(randn(3,1), rand(1)))';

  %x = 10*pi/180*randn(9,1);
  x  = zeros(9,1);
  x(4:end) = randn(6,1);
  x(1:3) = 200*pi/180*randn(3,1);
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

  [r, pred, drdx] = observe_box(y, x, q, dt, d0, r0);
  
  dr_dx = zeros(12,9);
  for i = 1:9
    dr_dx(:,i) = ( observe_box(y, x+dx*I9(:,i), q, dt, d0, r0) - r ) / dx;
  end
  
  if ( norm(drdx - dr_dx) > tol )
    disp('Test 1. Failed')
    cat(2, drdx, dr_dx)
    keyboard
  else
    disp('Test 1. OK')
  end


  N = 100;
  phi = linspace(0,2*pi,N)';
  sp = sin(phi);
  cp = cos(phi);
  sp2 = sin(phi+pi);
  cp2 = cos(phi+pi);
  z = zeros(size(sp));
  Y = [cp sp z sp cp z cp2 sp2 z sp2 cp2 z];
  
  nmarkers = 4;
  p0 = reshape(Y(1,:), 3, nmarkers );
  d0 = mean(p0, 2);
  r0 = cat(1, p0 - repmat(d0, 1, nmarkers), \
	   zeros(1, nmarkers));
    
  qfr = [0;0;0;1];
  
  [r, pred, drdx] = observe_box(Y(1,:)', zeros(9,1), qfr, dt, d0, r0);
  
  if norm(r) > tol
    disp('Test 2. Failed')
    r
    keyboard
  else
    disp('Test 2. OK')
  end
    
  [r, pred, drdx] = observe_box(Y(2,:)', zeros(9,1), qfr, dt, d0, r0);
  
  if norm(r) > tol
    disp('Test 2. Failed')
    r
    keyboard
  else
    disp('Test 2. OK')
  end
    
  

  keyboard
