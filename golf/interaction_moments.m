function [im, C] = interaction_moments(gm, states);
%  [im, C] = interaction_moments(gm, states);
% Function that calculates the interaction moment and Coriolis
% matrix for the kinematic model given the states sequence.
%
%% Input
%    gm       ->  model struct. Assumed to be single chain
%    states   ->  matrix (nsts x nfrs) with the estimated states of
%                 the model.
%% Output
%    im       <-  (nst x nfrs) matrix containing the interaction
%                 moments for each degree of freedom
%    C        <-  (nst x nst x nfrs) 3d array containing the
%                 Coriolis matrices


% Kjartan Halvorsen
% 2017-04-19


try

  % The endpoint path and mechanism Jacobian
  nst = size(states,1)/2;
  nfrs = size(states,2);

  [tws, g0, Mb] = flatten_km(gm);

  theta = states(1:nst, :);
  thetadot = states(nst+1:end, :);
  

  C = zeros(nst, nst, nfrs);
  im = zeros(nst, nfrs);
  for i=1:nfrs
    Ci = manipulator_coriolis_flattened(Mb, tws, g0, theta(:,i), ...
                                       thetadot(:,i));
    C(:,:,i) = Ci;
    im(:,i) = Ci*thetadot(:,i);
  end
  
catch
  keyboard
end
