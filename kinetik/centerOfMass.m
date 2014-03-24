function [com_all, com_segments] = centerOfMass(tws, coms, mass, states)
%  [com_all, com_segments] = com(tws, coms, mass, states)
%
% Computes the center of mass (CoM) of the coupled segments, as well as
% the CoM for each of the segments separately.
%
% Input
%    tws     -> nested cell array with twists
%    p0s     -> nested cell array with CoM for each segment, 
%               specified for the reference position.
%    mass    -> the masses for each of the segments.
%    states  -> the state vectors. (n x N) matrix, n is the number
%               of degrees of freedom. 
% Output
%    com_all <- The trajectory of the center of mass for the
%               set of segments. (3 x N) matrix.
%    com_segments <- The trajectory of the center of mass for each
%                    of the segments. (3 x nsgm x N)
  
% Kjartan Halvorsen
% 2003-07-29

totalmass = sum(mass);
massd = diag(mass);

[nst, nfr] = size(states);

ytmp = observe_mechanism_H( [states(:,1); zeros(nst,1)],1, ...
			    tws, coms);

nsgm = length(ytmp)/3;

com_all = zeros(3, nfr);
com_segments = zeros(3, nsgm, nfr);

for i=1:nfr
  y = observe_mechanism_H([states(:,i); zeros(nst,1)], 1, ...
			  tws, coms);
  comsegments = reshape(y, 3, nsgm);
  com_segments(:,:,i) = comsegments;
  
  comsegments = comsegments * massd;
 
  comt = sum(comsegments,2);
  com_all(:, i) = comt/totalmass;
end

  

  
