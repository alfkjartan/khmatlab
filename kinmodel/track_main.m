% Script for tracking a kinematic model.
% Test of new (faster?!) implementation.

% Kjartan Halvorsen
% 2011-10-26

% Global variables. Ugly, but works if careful.
global GL_xi0 GL_g GL_gs GL_p0 GL_p GL_pnames GL_pdep
%%global GL_xi0 GL_R GL_d0 GL_d GL_p0 GL_p GL_pnames GL_Rs GL_ds

debug = 1;

filterbandwidth = 1;  % Scales the process noise covariance
                      % matrix. Increase the value for higher
                      % bandwidth (less smoothing effect), lower
                      % for more smoothing.
skipframes = 4;   % The number of frames to skip in a downsampling of
                  % the data. 1 means no skipping.

datapth = ['.'];

trials = {'piano.tsv', 'piano_check'
          'accuracy.tsv', 'accuracy_check'};

%usetrials = (13:14);
usetrials = (2);

% Close all open plot windows
close all

% Names of markers to plot must match exactly with the names in the
% c3d file.
markers2plot = {'tip1'
	       'tip2'
	       'tip3'
	       'tip4'
	       'tip5'};
for i=1:length(markers2plot)
  markers2plot{i,2} = figure('Name', markers2plot{i,1});
end


% Names of joint angles to plot must match exactly with the names in the
% build_hand_model.m file.
angles2plot = {'mcp1 f/e'};
for i=1:length(angles2plot)
  angles2plot{i,2} = figure('Name', angles2plot{i,1});
end

% Now process the data
refdata = openmocapfile('', ...
			fullfile(datapth,'refnew.tsv'));

% Change from mm to m as unit
refdata{2} = refdata{2}/1000.0;

% Create model.
[hm] = build_hand_model(refdata);
%[hm] = build_simple_model(refdata);
segments = kinmodel2globals(hm); % Will set global variables and create
                                % set of segments

nstates = size(GL_xi0,3);
x0 = zeros(nstates, 1);
yt = ones(size(GL_p0));
[initnames, p0, p0vec] = prepare_mdata(hm.p0);

time_start = now();
for j=1:500
  [y] = observe_mechanism_H(x0, yt, hm.twists, p0);
end
oldtime = now() - time_start;

time_start = now();
for j=1:500
  [y] = observe_quick(x0, yt, segments);
end
newtime = now() - time_start;

disp(['timing resultat:'])
newtime
oldtime

disp(['Speed increase: ', sprintf("%f", oldtime/newtime)])


%save("-mat", "handmodel.mat",  "hm")

