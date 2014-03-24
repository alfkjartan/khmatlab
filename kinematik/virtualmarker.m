function vm = virtualmarker(refdata,mdata,segmarkers, ...
			    pointmarkers, nulldist)
%  vm = virtualmarker(refdata,mdata,segmarkers,pointmarkers)
% Calculates the position of the virtual marker indicated by the
% two pointer markers.
% The virtual marker is defined to lie on the line passing through
% the two pointer markers. The virtual marker lie at a distance
% from the first pointer marker in the direction away from the
% second pointer marker. The distance is given by the "null
% distance" (nulldist) minus the distance between the
% pointer markers.
% Input
%    refdata      ->   Marker data with reference position of
%                      segment.
%    mdata        ->   Marker data containing pointer markers
%    segmarkers   ->   Cell array with names of markers on the
%                      segment 
%    pointmarkers ->   Cell array with the two names of the markers
%                      on the pointer device.
%    nulldist     ->   distance between pointer markers defining
%                      zero depth.
% Output
%    vm           <-   The position of the virtual marker. For the
%                      segment markers in the position in refdata.

% Kjartan Halvorsen
% 2004-01-26

segmref = extractmeanmarkers(refdata, segmarkers);
segmm = extractmeanmarkers(mdata, segmarkers);

% Compute the rigid transformation taking the segment markers back
% to the reference position.
g_trf = soder(cat(1, segmm(:)', segmref(:)'));

pointm = extractmeanmarkers(mdata, pointmarkers);

pointvec = pointm(:,1) - pointm(:,2);
pdist = norm(pointvec);

vm = pointm(:,1) + (nulldist - pdist)* pointvec / pdist;

vm = g_trf * cat(1, vm, 1);
vm = vm(1:3);