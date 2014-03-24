function [Mn, Md, w, x] = euler_bernoulli_dynamics(F, EI, dt, N, K, a, b, c, m)

if (nargin == 0)
  do_unit_test();
else
  L = (2*a+b);
  my = m/L;
  dx = L/K;

  eyeK = eye(K);
  D = eyeK*6;
  D = D + diag(repmat(-4, K-1, 1), 1) + diag(repmat(-4, K-1, 1), -1) ...
      + diag(repmat(1, K-2, 1), 2) + diag(repmat(1, K-2, 1), -2);

  %% Fix finite difference approximations at the ends
  D(1,1:5) = [1 -4 6 -4 1];
  D(2,1:5) = [1 -4 6 -4 1];
  D(K-1,K-4:K) = [1 -4 6 -4 1];
  D(K,K-4:K) = [1 -4 6 -4 1];
  A12 = -EI/my/dx^4*D + 2/dt^2*eyeK;

  if (norm(A12) > 1)
    warning('Unstable difference equation')
    keyboard
  end

				%A = eye(2*K);
				%A(1:K, 1:K) = eyeK/dt^2;
				%A(1:K, K+1:2*K) = A12;
  
  %% Boundary conditions
  x = linspace(0, L, K);
  inda = round(a/L*K);
  indab = round((a+b)/L*K);

  %%A = cat(1, A, zeros(2, 2*K));
  %%A(2*K+1, inda) = 1;
  %%A(2*K+2, indaab) = 1;
  %%Ainv = inv(A);

  A = eye(K);
  A = cat(1, A, zeros(2, K));
  A(K+1, inda) = 1;
  A(K+2, indab) = 1;
  Ainv = inv(A'*A)*A';

  %% Right hand side. q/my. Step force
  q = zeros(K, 1);
  q(1) = F/my;
  rhs = zeros(K+2,1);

  wn = zeros(K,1);
  wnn = zeros(K, 1);

  w = zeros(K,N);

  for i=1:N
    wntmp = wnn;
    rhs(1:K) = dt^2*A12*wn - wnn + dt^2*q;

    wnn = wn;
    wn = Ainv*rhs;
    %%wn = A\rhs;
    
    w(:,i) = wnn;
  end

  Mn = zeros(N,1);
  Md = zeros(N,1);

end % if (nargin == 0)

function do_unit_test()
  a = 0.5;
  b = 0.6;
  c = 0.075;
  m = 1.1;

  F = 100;
  
  r = 15e-3; % radius of the shaft
  t = 1e-3; % thickness of wall
  II = pi*r^3*t; %% Thin tube approximation to second area of moment
  E = 200e9; % in GPa = 10^9 N / m^2
  EI = E*II;

  %% For stability, we need 
  dt = 1e-8;
  N = 100;
  K = 10;

  [Mn, Md, w, x] = euler_bernoulli_dynamics(F, EI, dt, N, K, a,  b, c, m);

  figure(1)
  clf

  grays = linspace(0, 0.7, N);
  for i = 1:10:N
    plot(x, w(:,i), 'color', repmat(grays(i), 1, 3));
    hold on
  end


  keyboard

