function W = mobility(km, states)
%%  W = mobility(km, states)
%% Returns the mobility tensor for the manipulator described by the model struct km.
%%
%% Input
%%    km         ->  The model struct
%%    states     ->  The sequence of states, (nsts x nfr)
%% Output
%%    W          <-  The mobility tensor (6 x 6 x nfr)

%% Kjartan Halvorsen
%% 2013-06-04


if nargin == 0
   do_unit_test();
else
  II = generalized_manipulator_inertia(km, states); % (nsts x nsts x nfr)

  
  [nsts, nfrs] = size(states);
  W = zeros(6, 6, nfrs);
  for i=1:nfrs
    state = states(:,i);
    Jse = spatial_manipulator_jacobian(km.twists, state); % (6 x nsts)
    W(:,:,i) = Jse * inv(II(:,:,i)) * Jse';
  end
end


function do_unit_test()

%% Construct model of the scara robot
l1 = 1;
l2 = 2;
m1 = 1;
m2 = 2;
m3 = 0.6;
m4 = 0.3;

sm = scara_robot_model(l1, l2, m1, m2, m3, m4);

states = [pi/4; pi/5; pi/6; pi/7];

Mtrue = scara_robot_inertia(sm, states(2));
Jstrue = spatial_manipulator_jacobian(sm.twists, states);

Wtrue = Jstrue * inv(Mtrue) * Jstrue'
W = mobility(sm, states)

