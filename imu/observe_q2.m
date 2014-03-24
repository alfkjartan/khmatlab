function [y,JJ]=observe_q2(x)
%  [y,JJ]=observe_q2(x)
% Measurement function for two imus fixed to the same object
%
% Input
%    x      ->   the current state. (n x N) vector.  x=[q12,q,w]

% Output
%    y     <-   the measurements [w1,w2]
%    JJ     <-   the Jacobian dy/dx

% Kjartan Halvorsen
% 2012-03-27
%
%

if (nargin == 0)
  do_unit_test();
else
  N = size(x,2);

  y = zeros(6,N);

  %% The angular velocity are given in the body frame.
  y(1:3,:) = x(9:11,:);

  for i=1:N
    %%y(4:6,i) = qtransv(x(9:11,i), x(1:4,i));
    y(4:6,i) = myqtransv(x(9:11,i), x(1:4,i));
  end

  if (nargout==2)
    %% Compute the jacobian. Assuming x has only one column
    %dw1/dqf = 0;
    %dw1/dq  = 0;
    dw1dw = eye(3);%dw1/dw  = I
    %dw2dqf = dqf/dqf*w*qconj(qf) + qf*w*d/dqf qconj(qf)
    dqfdqf = eye(4);
    dconjqfdqf = diag([-1 -1 -1 1]);
    dw2dqf = zeros(3,4);
    dw2dw = zeros(3,3);
    wq = [x(9:11); 0];
    qf = x(1:4);
    qfc = qconj(qf);
    wqqfc = qmult(wq,qfc);
    qfwq = qmult(qf, wq);
    for i=1:4
      ddqfi = qmult(dqfdqf(i,:),wqqfc) \
	  + qmult(qfwq,dconjqfdqf(i,:));
      dw2dqf(:,i) = ddqfi(1:3)';

      ddw = qmult(qf,qmult(dqfdqf(i,:), qfc));
      if (i<4)
	dw2dw(:,i) = ddw(1:3)';
      end
    end
    JJ = zeros(6,11);
    JJ(1:3,9:11) = dw1dw;
    JJ(4:6,1:4) = dw2dqf;
    JJ(4:6,9:11) = dw2dw;
  end
end

function vv = myqtransv(v,q)
  vq = cat(1, v, 0);
  vvq = qmult(qmult(q,vq), qconj(q));
  vv = vvq(1:3);

function do_unit_test()
  tol = 1e-5;
  
  disp("Unit test of function observe_q2")
  q12 = quaternion(randn(3,1),rand(1)*pi/2)';
  %%q12 = quaternion([1;0;0],pi/2)';
  %%q12 = [0;0;0;1];
  q = quaternion(randn(3,1),rand(1)*pi/2)';
  w = randn(3,1);
  %%w = zeros(3,1);
  %%w = 100*ones(3,1);
  x = cat(1, q12,q,w);
  
  [y,JJ] = observe_q2(x);

  xx = eye(11);
  dw = 1e-6;
  dydw = zeros(6,3);
  for i=9:11
    dydw(:,i-8) = (observe_q2(x+xx(:,i)*dw) - y) / dw;
  end


  if ( norm(JJ(:,9:11) - dydw) > tol )
    disp('Test 1. Failed')
    cat(2, JJ(:,9:11), dydw)
  else
    disp('Test 1. OK')
  end

  dydqf = zeros(6,4);
  for i=1:4
    dydqf(:,i) = (observe_q2(x+xx(:,i)*dw) - y) / dw;
  end

  if ( norm(JJ(:,1:4) - dydqf) > tol )
    disp('Test 2. Failed')
    cat(2, JJ(:,1:4), dydqf)
  else
    disp('Test 2. OK')
  end

  %%keyboard
