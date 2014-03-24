function [nallmarkers] = rotatec3dmarkers(markers, trfm, pg, names)
%   [nallmarkers] = rotatec3dmarkers(markers, trf [, pg, names])
% Rotates all the markers by the applied transformation 
%
% Input
%   markers      ->  marker data. (nfrs x nnmarks x 3)
%   trf          ->  transformation matrix (4 x 4)
%   pg           ->  parameter group. Optional
%   names        ->  cell array of names. Optional
% Output
%   nallmarkers     <-  new set of markers

% Kjartan Halvorsen
% 2005-03-08

[nfrs, nmrks, tre] = size(markers);

if (nargin==2)
  nmarkers = trfm * cat(1, reshape(permute(markers, [3 2 1]), ...
				   3, nfrs*nmrks), ...
			ones(1, nfrs*nmrks));
  nmarkers(4,:) = [];
  nallmarkers = permute(reshape(nmarkers, 3, nmrks, nfrs), [3 2 1]);
else
  labels = getc3dparam(pg, 'POINT', 'LABELS');

  [slask, nameind] = intersect(labels.data, names);

  if isempty(slask)
    warning(['Markers: ', join(names, ', '), '  not found in marker' ...
    ' set'])
    nallmarkers = [];
    return
  end
  
  nallmarkers = markers;
  
  nallmarkers(:, nameind, :) = ...
      rotatec3dmarkers(markers(:, namind, :), trfm);
end