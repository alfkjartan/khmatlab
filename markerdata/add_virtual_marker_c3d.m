function [nrefpos] =  add_virtual_marker_c3d(c3din, c3dout, c3dref,...
					     segmmarkers, newmarker, nref)
% add_virtual_marker_c3d(c3din, c3dout, segmmarkers, refpos, ...
%				newmarker, nref)
%
% Will add marker data to the c3d data and write to a new c3d
% file. Works by assuming that the new virtual marker is fixed to a
% segment for which enough (at least 3) markers exists.
% to 
%
% Input
%    c3din            ->  name of c3d file to read data from
%    c3dout           ->  name of c3d file to write
%    c3dref           ->  name of c3d file containing reference
%                         positions of the markers
%    segmmarkers      ->  cell array of markers fixed to the same
%                         segment.
%    newmarker        ->  name of the new marker
%    nref             ->  Optional. reference position of the new
%                         marker. (3 x 1) vector. If not provided,
%                         then either the marker is present in the
%                         c3d ref file, or, the position must be defined
%                         interactively.

% Kjartan Halvorsen
% 2007-04-17

% Revisions
% 2008-04-16   Added the possibility to use the reference position
%              of the marker in the reference c3d file

[markers, vidfrate, analogs, analogfrate, evnt, ParameterGroup,...
 ci, residerr] =  readC3D(c3din);
md = openmocapfile('',c3din);

mdref = openmocapfile('',c3dref);

refpos = extractmeanmarkers(mdref, segmmarkers);

if (nargin < 6)
  
  nref = extractmeanmarkers(mdref, {newmarker});
  
  if sum(nref) == 0

    disp(['You need to edit nref to provide the reference ',...
	  'position of the virtual markrer'])
    keyboard
  end
  
end

nrefpos = nref;
nmd = addmarkers(md, nref, {newmarker}, segmmarkers, mdref);

[nfrs, m3] = size(md{2});

origmarkers = permute(reshape((md{2})', [3 m3/3 nfrs]), [3 2 1]);

nmd = extractmarkers(nmd, {newmarker});

[allmarkers, npg, nci, nresids] = ...
    setc3dmarkers(permute(nmd, [1 3 2]),...
		  origmarkers, ParameterGroup, {newmarker}, ci, residerr);

writeC3D(allmarkers, vidfrate, analogs, analogfrate, evnt, ...
	 npg, nci, nresids, c3dout); 



