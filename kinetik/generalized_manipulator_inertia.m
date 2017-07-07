function [M, Mtrue] = generalized_manipulator_inertia(km, states)
%  M = generalized_inertia(km, states)
% Returns the generalized inertia for the manipulator described by the model struct km.
% See eq (4.19) and eq (4.29) in Murray, Li, Sastry 
%
% Input
%    km         ->  The model struct
%    states     ->  The sequence of states, (nsts x nfr)
% Output
%    M          <-  The generalized manipulator inertia matrix (nsts x nsts x nfr). That is
%                   the inertia matrix in joint space.

%% Kjartan Halvorsen
% 2013-05-29

if nargin == 0
   [M, Mtrue] = do_unit_test();
else
  [nsts, nfrs] = size(states);

  M = zeros(nsts, nsts, nfrs);

  Mlinks = link_inertia(km); % Returns a (6 x 6 x nsegms) matrix consisting of the nested
                             % inertia matrices in km.inertia. These must be given in the
                             % local coordinate system of each segment.
  nsegms = size(Mlinks,3);

  for i=1:nfrs
    state = states(:,i);
    Jsb = link_jacobians(km.twists, km.g0, state); %% Returns a (6 x nsts x nsgms) matrix
    for s=1:nsegms
      Jsbs = Jsb(:,:,s);
      M(:,:,i) = M(:,:,i) + Jsbs'*Mlinks(:,:,s)*Jsbs;
    end
  end

  %%keyboard
end

function [M, Mtrue] = do_unit_test()

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

%keyboard

M = generalized_manipulator_inertia(sm, states);

tol = 1e-12;

if norm(M - Mtrue) > tol
   disp('Test 1 failed')
   disp('Expected'), Mtrue
   disp('Found'), M
else
   disp('Test 1 OK')
end
