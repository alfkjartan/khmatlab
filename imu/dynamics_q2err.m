function [xx,JJ] = dynamics_q2err(x, q1, q12, dt)
%% [xx,JJ] = dynamics_q2err(x, q1, q12, dt)
%% Discrete time dynamics of single rigid body, orientation part. Two imus connected
%% to the same body
%% The state contains the orientation increment. The actual
%% (non-incremented) orientations are considered parameters to the model.
%% The state vector contains
%%    w12        <-   Increment for quaternion representing the rotation of imu2 wrt imu1
%%    w          <-   Angular velocity of RB. 
%% Output
%%    xx         <- predicted state 
%%    JJ         <- the jacobian dxx/dx. Works only for first column in x
  

%% Kjartan Halvorsen
%% 2012-03-29. 
%%

if nargin==0
  do_unit_test;
else
  xx = x;
  xx(1:3) = 0;
zeros(6,1);
  xx(4:5) = q = x(5:8,:);
  w = x(9:11,:);
  %aa = x(8:10,:);

  %% the orientation. No need to renormalize
  xx(5:8,:) = qpropagate(q,w,dt);

  %% the angular velocity. Random walk.
  %%  xx(5:7,:) = w + dt*aa; 

  if ( nargin == 2 )
    %% Calculate jacobian
    [slask, dqqdq, dqqdw] = qpropagate(q,w,dt);
  
    JJ = eye(11);
    JJ(5:8,5:8) = dqqdq;
    JJ(5:8,9:11) = dqqdw;
  end
end


function do_unit_test()
   tol = 1e-5;
  
  disp("Unit test of function dynamics_q2")
  q12 = quaternion(randn(3,1),rand(1)*pi/2)';
  %%q12 = quaternion([1;0;0],pi/2)';
  %%q12 = [0;0;0;1];
  q = quaternion(randn(3,1),rand(1)*pi/2)';
  w = randn(3,1);
  %%w = zeros(3,1);
  %%w = 100*ones(3,1);
  x = cat(1, q12,q,w);
  
  dt = 1;
  [x2,JJ] = dynamics_q2(x,dt);

  xx = eye(11);
  dw = 1e-6;
  dxxdw = zeros(11,3);
  for i=9:11
    dxxdw(:,i-8) = (dynamics_q2(x+xx(:,i)*dw,dt) - x2) / dw;
  end

 
  if ( norm(JJ(:,9:11) - dxxdw) > tol )
    disp('Test 1. Failed')
    cat(2, JJ(:,9:11), dxxdw)
  else
    disp('Test 1. OK')
  end

 
end  
