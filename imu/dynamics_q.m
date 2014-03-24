function [xx] = dynamics_q(x, dt)
% [xx] = dynamics_rbb(x, dt)
%% Like dynamics_rb, but using description in the body frame
% Discrete time dynamics of a rigid body using quaternial represenation
% of rotation.
% Should work for matrix valued x, with state vectors in columns.
% The state vector contains
%    q          <-   quaternion representing the rotation of the RB.
%                    Rotates back to the static frame
%    w          <-   Angular velocity of RB. 
%    alpha      <-   Angular accelerastion of RB.

% Kjartan Halvorsen
% 2012-03-29. 
%

if nargin==0
  do_unit_test;
else
  xx = x;
  q = x(1:4,:);
  w = x(5:7,:);
  %aa = x(8:10,:);

  %% the orientation. No need to renormalize
  xx(1:4,:) = qpropagate(q,w,dt);

  %% the angular velocity. Random walk.
  %%  xx(5:7,:) = w + dt*aa; 
end


function do_unit_test()
  
end  
