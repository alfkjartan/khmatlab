function g = get3dofMotion( mdata, refdata, proxmarkers, distmarkers, ...
				   jc)
%  
%  g = get3dofMotion( mdata, refdata, proxmarkers, distmarkers, ...
%				   jc)
% Computes the motion of the distal segment, given the motion of
% the segment proximal, and a 3dof ball and socket joint with known
% position.
% Input
%    mdata         -     Marker data. (2 x 1) cell array: {attr,
%                        md}
%    refdata       -     Marker data. (2 x 1) cell array: {attr,
%                        md}. Reference positon of markers.
%    proxmarkers   -     Cell array of marker names for the
%                        proximal segment.
%    distmarkers   -     Cell array of marker names for the distal segment.
%    jc            -     position of the 3dof joint center.
% Output
%    g             -     (4 x 4 x nfr) transformation matrices. 

% Kjartan Halvorsen
% 2003-11-17

if (length(distmarkers) < 2)
  error('Not enough markers to compute motion')
end

proxm = extractmarkers(mdata, proxmarkers);
proxmref = extractmeanmarkers(refdata, proxmarkers);

distm = extractmarkers(mdata, distmarkers);
distmref = extractmeanmarkers(refdata, distmarkers);

gprox = getMotion(cat(1, proxmref(:)', proxm));

nfr = size(gprox, 1);
jcp = movePoints(gprox, jc);
jcp = reshape(jcp, 3, nfr);

% Add the joint center to the data for the distal marker.
% Then run the getMotion function
g = getMotion(cat(2, cat(1, distmref(:)', distm), jcp'));

g=g(2:end, :);


  