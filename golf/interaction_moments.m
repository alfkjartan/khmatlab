function [im, C, dofnames] = interaction_moments(gm, states, dofs2study);
%  [im, C] = interaction_moments(gm, states [, dofs2study]);
% Function that calculates the interaction moment and Coriolis
% matrix for the kinematic model given the states sequence.
%
%% Input
%    gm       ->  model struct. Assumed to be single chain
%    states   ->  matrix (nsts x nfrs) with the estimated states of
%                 the model.
%    dofs2study -> list of names of dofs for which to calculate
%                  interaction moments (optional) 
%% Output
%    im       <-  (nst x nfrs) matrix containing the interaction
%                 moments for each degree of freedom to study. The sign
%                 corresponds to the sign of a muscle torque acting over
%                 the same joint.
%    C        <-  (nst x nst x nfrs) 3d array containing the
%                 Coriolis matrices. These have sign corresponding to the
%                 coriolis term in the left hand side of the equations of
%                 motion M(q)\ddot{q} + C(q, \dot{q})\dot{q} + G(q) = tau


% Kjartan Halvorsen
% 2017-04-19


%try

  % The endpoint path and mechanism Jacobian
  nst = size(states,1)/2;
  nfrs = size(states,2);

  [tws, g0, Mb, gcnames] = flatten_km(gm);
  
  if nargin < 3
      dofs2study = gcnames(:,1);
  end
  
  [dofnames_unsorted, thetas2study_unsorted] = intersect(gcnames(:,1), dofs2study);

  [thetas2study, ids] = sort(thetas2study_unsorted);
  dofnames = dofnames_unsorted(ids);


  theta = states(1:nst, :);
  thetadot = states(nst+1:end, :);
  
  ndofs2study = length(thetas2study);
  
  C = zeros(ndofs2study, ndofs2study, nfrs);
  im = zeros(ndofs2study, nfrs);
  for i=1:nfrs
    Ci = manipulator_coriolis_flattened(Mb, tws, g0, theta(:,i), ...
                                       thetadot(:,i), []);
    C(:,:,i) = Ci(thetas2study, thetas2study);
    im(:,i) = Ci(thetas2study, :)*thetadot(:,i);
  end
  
  im = -im;
%catch
%  keyboard
%end
