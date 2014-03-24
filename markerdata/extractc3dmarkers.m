function m = extractc3dmarkers(markers, pg, names)
% m = extractc3dmarkers(markers, pg, names)
% Extracts the trajectories for the marker names provided.
%
% Input
%   markers         ->  marker data as returned by readC3D (nfrs x
%                       nmarkers x 3)
%   pg              ->  parameter group
%   names           ->  cell array of names.
% Output
%   m               <-  the marker trajectories 
%                       (nfrs x length(names) x 3)

% Kjartan Halvorsen
% 2005-03-01

labels = getc3dparam(pg, 'POINT', 'LABELS');

[slask, distind, nameind] = intersect(labels.data, names);
    
% Note that the indices are sorted. Fix this
[slask, sortind] = sort(nameind);
distind = distind(sortind);
    
[nfrs, nmarkers, tre] = size(markers);

m = zeros(nfrs, length(names), 3);
for i=1:length(distind)
  m(:,i,:) = markers(:,distind(i),:);
end


