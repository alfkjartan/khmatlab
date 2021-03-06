function [ac, ep_vel] = anglecontributions_test(gm, states);
% Function to test the anglecontribution function on a simpler
% model with only a single chain. 
%  ac = anglecontributions(gm, states, endpoint);
% Function that calculates each angles contribution to the end point
% velocity. This is taken as the velocity of the CoM of the
% last segment. See build_golf_model for the definition of
% this point.
% Let p(theta) = g(theta)p_0 be the function that gives the
% position of the end point as a function of the joint angles
% (and other generalized coordinates). The contribution is
% then calculated as the derivative of this function w.r.t
% each joint angle. The resulting expression expresses the
% infinitesimal change in the end point position for an
% infinitesimal change in the joint angle, i.e. a vector. This
% vector is multiplied with the actual joint angle velocity,
% and then projected onto the direction of the actual velocity
% of the end point.
%
% Input
%    gm       ->  model struct
%    states   ->  matrix (nsts x nfrs) with the estimated states of
%                 the model.
% Output
%    ac       ->  (nst x nfrs) matrix containing the contributions
%                 of each joint angle ( generalized coordinate).

% Kjartan Halvorsen
% 2009-06-29


try
  % The endpoint path and mechanism Jacobian
  nst = size(states,1)/2;
  % Assign the CoM initial points to the models initial points p0,
  % to obtain the Jacobian for these.
  gm.p0 = gm.CoM;
  [ep, epnames, slask1, slask2,J] = sim_model(gm,states(1:nst,:), 'CoM');

  %keyboard
  
  % Hard-coded sampling frequency of 120Hz
  ep_vel = centraldiff(ep, 120);

  nfrs = size(states,2);

  ac = zeros(nst,nfrs);
  for i=1:nfrs
    ep_vel_dir1 = ep_vel(i,1:3) / norm(ep_vel(i,1:3));
    %ep_vel_dir1 = [0 1 0];
    %ep_vel_dir2 = [0 1 0];

    Jscaled = J(:,:,i)*diag(states(nst+1:end,i));
    %Jscaled = J(:,:,i);
    ac(:,i) = (ep_vel_dir1*Jscaled(1:3,:))';
    
  end

catch
  keyboard
end
