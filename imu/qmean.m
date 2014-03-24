function [qt, E, its] = qmean(Q,W, thr)
%% [qm,E] = qmean(Q,W)
%% Computes the mean of a set of quaternions - the columns of Q.

%% Kjartan Halvorsen
%% 2012-03-20

if (nargin == 0)
  do_unit_test();
else

  MAXITS = 10;
  if nargin < 3
    thr = 1e-4; % Break iteration when change is less than tolerance
  end

  if nargin < 2
    W = [];
  end

  tol = 1e-10; % Rotations below this are identity rotations

  %% starting guess
  qt = Q(:,1);

  N = size(Q,2);
  E = zeros(3,N);

  adjust = 1;
  its = 0;
  while (adjust > thr)
    its = its + 1;
    if (its > MAXITS)
      its = -1;
      return
    end
    if (its>1)
      etn = et / adjust;
      qt = qmult(cat(1, etn*sin(adjust/2), cos(adjust/2)), qt)';
      %%qt = qmult(qt,cat(1, etn*sin(adjust/2), cos(adjust/2)))';
    end
    qtinv = qinv(qt);
    for i=1:N
      ei = qmult(qtinv,Q(:,i));
      %%ei = qmult(Q(:,i),qtinv);
      th = 2 * acos(ei(4));
      if (abs(th) < tol)
	E(:,i) = zeros(3,1);
      else
	E(:,i) = ei(1:3)/sin(th/2)*th;
      end
    end

    if (isempty(W) == 0)
      et = zeros(3,1);
      for i=1:N
	et = et + W(i)*E(:,i);
      end
    else
      et = mean(E,2);
    end
    adjust = norm(et);
  end

end

function do_unit_test()
  warning('error', 'Octave:divide-by-zero')
  disp('Unit test of function qmean ')

  thr = 1e-5
  tol = 1e-4;

  n = 30;
  th = linspace(-pi/3,pi/3,n);

  w = randn(3,1);
  w = w/norm(w);
  
  Q = zeros(4,n);

  for i=1:n
    Q(:,i) = quaternion(w,th(i));
  end

  W = repmat(1/n, n,1);

  [qm1,E1,its1] = qmean(Q, [], thr);
  [qm2,E2,its2]  = qmean(Q,W, thr);

  disp(sprintf("Number of iterations: %d and %d", its1, its2))

  if (norm(qm1-qm2) > tol)
    disp('Test1: Failed')
    disp('Unexpected result:')
    disp('qm1='), disp(qm1)
    disp('qm2='), disp(qm2)
  else
    disp('Test1: OK')
  end

  if (abs(qm1(4) - 1) > tol)
    disp('Test2: Failed')
    disp('Unexpected result:')
    disp('Expected identity quaternion.')
    disp('Found '), disp(qm1)
  else
    disp('Test2: OK')
  end

  if (abs(norm(qm1) - 1) > tol)
    disp('Test3: Failed')
    disp('Unexpected result:')
    disp('Expected unit quaternion.')
    disp('Found '), disp(qm1)
  else
    disp('Test3: OK')
  end

  P=eye(3)*pi/6;
  w1 = randn(3,1);
  th1 = pi/3;
  q1 = quaternion(w1,th1)';
  q2 = [0;0;0;1];
  [X1,Wm,Wc] = sigmaq(P, q1, 1);
  [X2,Wm,Wc] = sigmaq(P, q2, 1);

  [qm1,E1,its1] = qmean(X1, Wm);
  [qm2,E2,its2] = qmean(X2, Wm);

  disp(sprintf("Test4: Iterations=%d", its1));
  if (norm(qm1 - q1) > tol)
    disp('Test4: Failed')
    disp('Unexpected result:')
    disp('Expected '), disp(q1)
    disp('Found '), disp(qm1)
  else
    disp('Test4: OK')
  end
 
  disp(sprintf("Test5: Iterations=%d", its2));
  if (norm(qm2 - q2) > tol)
    disp('Test5: Failed')
    disp('Unexpected result:')
    disp('Expected '), disp(q2)
    disp('Found '), disp(qm2)
  else
    disp('Test5: OK')
  end
 

  %% Test non-diagonal covariance matrix
  S = randn(3,3);
  P = S'*S;
  w = randn(3,1);
  w = w / norm(w);
  th = pi/2.5;
  q0 = quaternion(w,th)';
  [X,Wm,Wc] = sigmaq(P,q0, 1);
  
  [qm,EE,its] = qmean(X, Wm);

  disp(sprintf("Test6: Iterations=%d", its));
  if (norm(qm - q0) > tol)
    disp('Test6: Failed')
    disp('Unexpected result:')
    disp('Expected '), disp(q0)
    disp('Found '), disp(qm)
  else
    disp('Test6: OK')
  end
 
  Pe = zeros(3,3);
  for i=1:size(EE,2)
    Pe = Pe + Wc(i)*EE(:,i)*EE(:,i)';
  end

  if (norm(Pe - P) > tol)
    disp('Test7: Failed')
    disp('Unexpected result:')
    disp('Expected '), disp(P)
    disp('Found '), disp(Pe)
  else
    disp('Test7: OK')
  end
 
  

  