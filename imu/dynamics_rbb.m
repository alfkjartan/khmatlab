function [xx] = dynamics_rbb(x, dt)
% [xx] = dynamics_rbb(x, dt)
%% Like dynamics_rb, but using description in the body frame
% Discrete time dynamics of a rigid body using quaternial represenation
% of rotation.
% Should work for matrix valued x, with state vectors in columns.
% The state vector contains
%    d          <-   the estimated translation of the rigid body. In the
%                    static body frame
%    vel        <-   velocity of the RB
%    acc        <-   acceleration of the RB
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

  n = size(x,1);

  if (n==19)
    d = x(1:3,:);
    v = x(4:6,:);
    a = x(7:9,:);
    q = x(10:13,:);
    w = x(14:16,:);
    aa = x(17:19,:);

    %% the acceleration. No dynamics, it's a random walk
    %%xx(7:9,:) = xx(7:9,:);

    %% the orientation. No need to renormalize
    q2 = qpropagate(q,w,dt);    
    xx(10:13,:) = q2;

    %% the velocity. Rotate acceleration back to static frame. Use both
    %% last and current orientation estimate
    invert = 1;
    %a1 = qrot(a,q, invert);
    %a2 = qrot(a,q2, invert);
    %v2 = v + 0.5*dt*a1 + 0.5*dt*a2;
    v2 = v + dt*a;
    xx(4:6,:) = v2;

    %% the position
    %xx(1:3,:) = d + 0.5*dt*v + 0.5*dt*v2 + 0.25*dt^2*a1 + 0.25*dt^2*a2;
    xx(1:3,:) = d + dt*v + 0.5*dt^2*a;

    %% the angular velocity
    xx(14:16,:) = w + dt*aa;

    %% the angular acceleration 
    %%xx(17:19,:) = xx(17:19,:);
  elseif (n==13)
    d = x(1:3,:);
    v = x(4:6,:);
    q = x(7:10,:);
    w = x(11:13,:);

    %% the position
    xx(1:3,:) = d + dt*v;

    %% the orientation. No need to renormalize
    xx(7:10,:) = qpropagate(q,w,dt);
  else
    error(sprintf("Strange size of state vector: %d", n))
  end
end

function do_unit_test()
  
end  
