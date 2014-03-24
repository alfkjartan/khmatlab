function [ac, ac_l, ac_r, ep_vel_mean, ep_vel_l, ep_vel_r] ...
    = anglecontributions(gm, states);
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

% Revisions
% 2010-03-21  Added output ac_l and ac_r, which are the
%             contributions for the left and right sides separated,
%             from root segment to end point

try
  % The endpoint path and mechanism Jacobian
  nst = size(states,1)/2;
  % Assign the CoM initial points to the models initial points p0,
  % to obtain the Jacobian for these.
  gm.p0 = gm.CoM;
  [ep, epnames, slask1, slask2,J] = sim_model(gm,states(1:nst,:), 'CoM');
  % The above function will return two end points, due to the fact
  % that the club is at the end of two chains. The left is the first.

  % Hard-coded sampling frequency of 120Hz
  ep_vel = centraldiff(ep, 120);

  ep_vel_l = ep_vel(:,1:3);
  ep_vel_r = ep_vel(:,4:6);
  
  %keyboard
  
  % The velocities are almost the same for the left and right end
  % point. Let's therefore use the mean velocity as the assumed
  % true endpoint velocity.  For a different approach, registering
  % the max contrib to either left or right hand side endpoint look
  % at the backup file anglecontributions.m.100321
  
  ep_vel_mean = cat(2, ...
		    mean(ep_vel(:,[1 4]), 2), ...
		    mean(ep_vel(:,[2 5]), 2), ...
		    mean(ep_vel(:,[3 6]), 2));
  
  nfrs = size(states,2);

  ac_l = zeros(nst,nfrs);
  ac_r = zeros(nst,nfrs);
  for i=1:nfrs
    ep_vel_dir = ep_vel_mean(i,1:3) / norm(ep_vel_mean(i,1:3));
    ep_vel_dir1 = ep_vel(i,1:3) / norm(ep_vel(i,1:3));
    ep_vel_dir2 = ep_vel(i,4:6) / norm(ep_vel(i,4:6));

    Jscaled = J(:,:,i)*diag(states(nst+1:end,i));

    %ac_l(:,i) = (ep_vel_dir*Jscaled(1:3,:))';
    %ac_r(:,i) = (ep_vel_dir*Jscaled(4:6,:))';
    ac_l(:,i) = (ep_vel_dir1*Jscaled(1:3,:))';
    ac_r(:,i) = (ep_vel_dir2*Jscaled(4:6,:))';
  end
 
  ac = ac_l;
  rightside_gc = find(sum(ac_l,2) == 0);
  ac(rightside_gc,:) = ac_r(rightside_gc,:);
  
catch
  keyboard
end
