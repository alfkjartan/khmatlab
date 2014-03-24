%% Script for analyzing validation data from 2012-02-07

%% Kjartan Halvorsen
%% 2012-02-21

datapth = '/home/kjartan/Dropbox/projekt/golf/data/validering120207';

floord = openmocapfile('', fullfile(datapth,'floor.tsv'));
floormarker = extractmeanmarkers(floord, '?');

pointerd  = openmocapfile('', fullfile(datapth,'pointer.tsv'));
pointernames = {'pointer_tip'; 'pointer_handle'};
ptip = extractmeanmarkers(pointerd, pointernames{1});
phandle = extractmeanmarkers(pointerd, pointernames{2});
v = floormarker - ptip;
%%pointerdist = norm(v) - 5e-3;  % subtract radius of marker
%%pointerdist = norm(v) + 5e-3;  %  add radius of marker
pointerdist = norm(v) + 8e-3;  %  add radius of marker
%% Note that the floormarker was placed on the floor and measured separately. After that
%% the marker was removed and the pointer was positioned to point at the point 
%% where the marker was.

ref = fullfile(datapth, 'ref.tsv');
pliancepos = cell(12,1);
pliancepos{1} = fullfile(datapth, 'pos31.tsv');
pliancepos{2} = fullfile(datapth, 'pos32.tsv');
pliancepos{3} = fullfile(datapth, 'pos33.tsv');
pliancepos{4} = fullfile(datapth, 'pos34.tsv');
pliancepos{5} = fullfile(datapth, 'pos91.tsv');
pliancepos{6} = fullfile(datapth, 'pos92.tsv');
pliancepos{7} = fullfile(datapth, 'pos93.tsv');
pliancepos{8} = fullfile(datapth, 'pos94.tsv');
pliancepos{9} = fullfile(datapth, 'pos131.tsv');
pliancepos{10} = fullfile(datapth, 'pos132.tsv');
pliancepos{11} = fullfile(datapth, 'pos133.tsv');
pliancepos{12} = fullfile(datapth, 'pos134.tsv');

mpdata = {fullfile(datapth, 'weight1000_1.tsv'),...
	  fullfile(datapth, 'Tinmark_Fredrik_3.asc')};
%%mpdata = {fullfile(datapth, 'weight1000_1.tsv'),...
%%	  './test.asc'};
get_hand_force_static(ref, pliancepos, mpdata, pointerdist)
