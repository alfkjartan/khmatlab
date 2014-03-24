function md = qtm_export2md_format(qtmm) 
%  md = qtm_export2md_format(qtmm) 
%
% Transforms the exported marker data from qtm's format when
% exporting directly to matlab to the format used in most of my
% code:
% md = { attributes, data }, with 
% attributes = {key, value; key, value;...}, and
% data = nfr x 3*nmarkers
%
% OBS. Does not currently support analog data .

% Kjartan Halvorsen
% 2008-01-23

nfrs = size(qtmm.Trajectories.Labeled.Data, 3);
nmrkrs = qtmm.Trajectories.Labeled.Count;

attributes = { 'NO_OF_FRAMES', num2str(nfrs) ;...
	       'NO_OF_MARKERS', num2str(nmrkrs); ...
	       'FREQUENCY', num2str(qtmm.FrameRate);...
	       'MARKER_NAMES', qtmm.Trajectories.Labeled.Labels' };

%keyboard

dta = reshape(permute(qtmm.Trajectories.Labeled.Data, [3 2 1]), ...
	      [nfrs 3*nmrkrs]);

md = {attributes, dta};

