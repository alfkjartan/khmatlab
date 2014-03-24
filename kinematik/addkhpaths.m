function addkhpaths(extras)
% addkhpaths(extras)
% extras is an array of strings indicating extra directories to include
% in path

% Kjartan Halvorsen
% 2011-06-27

pth = fileparts(mfilename("fullpath"));

addpath(fullfile(pth,"util"));
addpath(fullfile(pth,"ekf"));
addpath(fullfile(pth,"kinematik"));
addpath(fullfile(pth,"kinetik"));
addpath(fullfile(pth,"markerdata"));
addpath(fullfile(pth,"signal"));
%addpath(fullfile(pth,"thirdpart"));

if (nargin>0)
for i=(1:length(extras))
  addpath(fullfile(pth, extras{i}));
endfor
endif

        
	
