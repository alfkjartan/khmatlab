function [xx] = dynamics_rb(x, dt)
% [xx] = dynamics_rb(x, dt)
% Discrete time dynamics of a rigid body using quaternial represenation
% of rotation.
% Should work for matrix valued x, with state vectors in columns.
% The state vector contains
%    d          <-   the estimated translation of the rigid body. The
%                    local coordinate of the RB coincides with the
%                    static lab frame for the first frame of data.
%    vel        <-   velocity of the RB, expressed in the static frame.
%    acc        <-   acceleration of the RB, expressed in the static frame.
%    q          <-   quaternion representing the rotation of the RB.
%                    Rotates a point in the first frame to subsequent positions.
%    w          <-   Angular velocity of RB. Expressed in the static
%                    frame.
%    e          <-   noise vector.
%    alpha      <-   Angular accelerastion of RB. Expressed in the
%                   static frame.

% Kjartan Halvorsen
% 2012-03-20. 
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

    %% the position
    xx(1:3,:) = d + dt*v + 0.5*dt^2*a;

    %% the velocity
    xx(4:6,:) = v + dt*a;

    %% the acceleration 
    %%xx(7:9,:) = xx(7:9,:);

    %% the orientation. No need to renormalize
    xx(10:13,:) = qpropagate(q,w,dt);

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
