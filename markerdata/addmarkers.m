function nmd = addmarkers(mdata, npoints, nmarkernames, markernames, refdata)
%  nmd = addmarkers(mdata, npoints, nmarkernames, markernames, refdata)
%
% Adds (virtual) points to the marker data. The points move
% rigidly with the markers given in markernames. The position of
% the points is provided for the reference position given in
% refdata.
%
% Input
%    mdata         ->  Marker data {attr, md}
%    npoints       ->  3D points [3 x np]
%    nmarkernames  ->  Cell array with marker names
%    markernames   ->  Cell array with marker names
%    refdata       ->  Marker data {attr, md} in reference position
% Output
%    nmd           <-  Marker data {attr, md}.

% Kjartan Halvorsen
% 2004-10-13

% Get the rigid body transformation
T = getMotion(refdata, mdata, markernames);

% Generate 3D paths
newp = movePoints(T, npoints);
np = size(newp,2);
nfr = size(newp,3);

% Add to markerdata
mnames = getvalue(mdata{1}, 'MARKER_NAMES');
mnames = cat(1, mnames, nmarkernames);
nmd = {putvalue(mdata{1}, 'MARKER_NAMES', mnames), ...
       cat(2, mdata{2}, (reshape(newp, 3*np, nfr))')};



