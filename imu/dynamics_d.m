function [xx, JJ] = dynamics_d(x, dt)
% [xx, JJ] = dynamics_d(x, dt)
%% Like dynamics_rb, but using description in the body frame
% Discrete time dynamics of a rigid body using quaternial represenation
% of rotation.
% Should work for matrix valued x, with state vectors in columns.
% The state vector contains
%    d          <-   the estimated translation of the rigid body. In the
%                    static body frame
%    vel        <-   velocity of the RB
%    acc        <-   acceleration of the RB
% Output
%%   xx         <- predicted state
%%   JJ         <- Jacobian dxx/dx

% Kjartan Halvorsen
% 2012-03-29. 
%

if nargin==0
  do_unit_test;
else
  xx = x;

  d = x(1:3,:);
  v = x(4:6,:);
  a = x(7:9,:);

  xx(1:3,:) = d + dt*v;

  xx(4:6,:) = v + dt*a;

  if (nargout > 1)
    Idt = dt*eye(3);
    JJ = eye(9);
    JJ(1:3, 4:6) = Idt;
    JJ(4:6, 7:9) = Idt;
  end
end

function do_unit_test()

