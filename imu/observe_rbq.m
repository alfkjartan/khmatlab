function [y]=observe_rbq(x, p0, d0)
%  [y]=observe_rb(x, p0)
% Measurement function for a rigid body with orientation represented by
% a quaternion in the first 4 elements of the state vector.
%
% Input
%    x      ->   the current state. (n x N) vector. 
%    p0      ->  [3 x m] matrix containing the initial position of the markers
% Output
%    y     <-   the current marker positions. (3*m x N) vector.

% Kjartan Halvorsen
% 2012-03-20
%
%

if (nargin == 0)
  do_unit_test();
else
  N = size(x,2);
  nmarkers = size(p0,2);

  y = zeros(nmarkers*3, N);

  for sp = 1:N
    d = x(5:7,sp);
    q = x(1:4, sp);

    %% The position of the markers
    for j = 1:nmarkers
      y((j-1)*3+1:j*3,sp) = d0 + d + qtransv(p0(:,j), q);
    end
  end
end

function do_unit_test()
  disp("Unit test for function observe_rb")

  tol = 1e-10;

  p0 = randn(24,1);
  p00 = reshape(p0, 3, 8);
  d0 = mean(p00, 2);
  p00 = p00 - repmat(d0, 1, 8);
 
  d = randn(3,1);

  p11 = p00 + repmat(d+d0, 1, 8);
  p1 = p11(:);

  p22 = p00;
  slask = p22(1,:);
  p22(1,:) = -p22(2,:);
  p22(2,:) = slask;
  p22 = p22 + repmat(d0, 1, 8);
  p2 = p22(:);

  x0 = zeros(19,1);
  qi = 1;
  x0(qi+3) = 1;

  x1 = x0;
  x1(5:7) = d;

  x2 = x0;
  x2(1:4) = quaternion([0;0;1], pi/2);

  y = observe_rbq(x0, p00, d0);

  if (norm(y - p0) > tol)
    disp(sprintf("Test 1: Failed. Norm = %2.8f", norm(y-p0)))
    disp("[Expected found] = "), disp(cat(2, y, p0))
  else
    disp('Test 1: OK')
  end


  y = observe_rbq(x1, p00, d0);
  if (norm(y - p1) > tol)
    disp(sprintf("Test 2: Failed. Norm = %2.8f", norm(y-p0)))
    disp("[Expected found] = "), disp(cat(2, y, p1))
  else
    disp('Test 2: OK')
  end

  y = observe_rbq(x2, p00, d0);
  if (norm(y - p2) > tol)
    disp(sprintf("Test 3: Failed. Norm = %2.8f", norm(y-p0)))
    disp("[Expected found] = "), disp(cat(2, y, p2))
  else
    disp('Test 3: OK')
  end

  
