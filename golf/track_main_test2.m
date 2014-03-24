% function track_main
% Builds, tracks and computes joint velocities contribution to
% the endpoint velocity.

% Kjartan Halvorsen
% 2009-06-26

filterbandwidth = 1;  % Scales the process noise covariance
                      % matrix. Increase the value for higher
                      % bandwidth (less smoothing effect), lower
                      % for more smoothing.

datapth = ['C:\Documents and Settings\Kjartan Halvorsen\My Documents',...
	   '\projekt\gih\golf\SAS'];
% Initials, folder name and reference mat-filee
fps = {'DP', 'DP', 'wedgemodel_dp.mat'
       'FP', 'FP', 'wedgemodel_fp.mat'
       'LJ', 'LJ', 'wedgemodel_lj.mat'
       'MH', 'MH', 'wedgemodel_mh.mat'
       'MW', 'MW', 'wedgemodel_mw.mat'
       'SS', 'SS', 'wedgemodel_ss.mat'};

trials = {'File251.c3d'
	  'File252.c3d'
	  'File253.c3d'
	  'File401.c3d'
	  'File402.c3d'
	  'File403.c3d'
	  'File551.c3d'
	  'File552.c3d'
	  'File553.c3d'
	  'File701.c3d'
	  'File702.c3d'
	  'File703.c3d'
	  'Filefull1.c3d'
	  'Filefull2.c3d'
	  'Filefull3.c3d'};

usefps = (1); % Select all or part of data to process
%usetrials = (13:14);
usetrials = (13);

% Close all open plot windows
close all

% Names of markers to plot must match exactly with the names in the
% c3d file.
markers2plot = {'L_Hand_1'
	       'R_Hand_1'};
for i=1:length(markers2plot)
  markers2plot{i,2} = figure('Name', markers2plot{i,1});
end


% Names of joint angles to plot must match exactly with the names in the
% build_golf_model.m file.
angles2plot = {'pelvis rotation'
	      'trunk rotation'};
for i=1:length(angles2plot)
  angles2plot{i,2} = figure('Name', angles2plot{i,1});
end

% Names of joint angles to plot must match exactly with the names in the
% build_golf_model.m file.
anglecontribs2plot = {'pelvis x'
		    'pelvis y'
		    'pelvis rotation'
		    'trunk rotation'};
for i=1:length(anglecontribs2plot)
  anglecontribs2plot{i,2} = figure('Name', ...
				   [anglecontribs2plot{i,1},' contrib']);
end


% Now process the data
for fp=usefps
  
  refdata = load(fullfile(datapth,fps{fp,2},fps{fp,3}));
  

  for tr=usetrials
    
    % Read the marker data
    filestr = fullfile(datapth,fps{fp,2},trials{tr});
    mdata = openmocapfile('', filestr);

    % Create model.
    % The trial data are needed because the position of the club
    % with respect to either hand is taken from the first frame
    % (address) of the trial file.
    [gmleft, gmright, epleftstr, eprightstr] = ...
	build_test2_model(refdata,mdata);
    
    % Edit below to switch analysis between left and right hand
    % side
    gm = gmright;
    epstr = eprightstr;
    gm = gmleft;
    epstr = epleftstr;
    
    [states, dataframes] = track_golf_model(gm, mdata, filterbandwidth);
    [states_l, dataframes_l] = track_golf_model(gmleft, mdata, filterbandwidth);
    [states_r, dataframes_r] = track_golf_model(gmright, mdata, filterbandwidth);
    nst = size(states,1)/2;
    
    % Simulate model, generate trajectories of joint centers.
    % plot markers and check the residuals
    %[msim,simnames, jcs] = sim_model(gm, states, 'jcs');
     [msim,simnames, jcs, jcnames] = sim_model(gm, states(1:nst,:), 'jcs');
     [msim,simnames, com, comnames] = sim_model(gm, states(1:nst,:), 'CoM');
     
     % Plot results
     y_observations=extractmarkers(mdata, simnames);
     y_observations(find(y_observations==0)) = NaN;

     
     mse = plotmarkers(y_observations(find(dataframes),:), ...
		       simnames,...
		       msim, simnames, markers2plot, {'data', 'model'});
    
     convert_radians = 1;
     plotangles(states, gm.gcnames(:,1), angles2plot, convert_radians);

     % Calculate each angles contribution to the end point
     % velocity. This is taken as the velocity of the CoM of the
     % club segment. See build_golf_model for the definition of
     % this point.
     % Let p(theta) = g(theta)p_0 be the function that gives the
     % position of the end point as a function of the joint angles
     % (and other generalized coordinates). The contribution is
     % then calculated as the derivative of this function w.r.t
     % each joint angle. The resulting expression expresses the
     % infinitesimal change in the end point position for an
     % infinitesimal change in the joint angle, i.e. a vector. This
     % vector is multiplied with the actual joint angle velocity,
     % and then projected onto the direction of the actual velocity
     % of the end point.
     [ang_contribs, model_endpoint_vel] = anglecontributions_test(gm, states);
     [ang_contribs_l, model_endpoint_vel_l] = anglecontributions_test(gmleft, states_l);
     [ang_contribs_r, model_endpoint_vel_r] = anglecontributions_test(gmright, states_r);
     %plotangles(ang_contribs, gm.gcnames(:,1), anglecontribs2plot, 0, 'm/s');

     % Debug angle contributions. Plot raw end velocity together
     % with sum of contributions and model end velocity.
     model_endpoint_vel_magn_l = sqrt(sum(model_endpoint_vel_l.^2, 2));
     model_endpoint_vel_magn_r = sqrt(sum(model_endpoint_vel_r.^2, 2));
     
     endpoint = extractmarkers(mdata, epstr);
     endpoint_vel = centraldiff(endpoint,120);
     endpoint_vel_magn = sqrt(sum(endpoint_vel.^2, 2));
     endpointleft = extractmarkers(mdata, epleftstr);
     endpoint_vel_l = centraldiff(endpointleft,120);
     endpoint_vel_magn_l = sqrt(sum(endpoint_vel_l.^2, 2));
     endpointright = extractmarkers(mdata, eprightstr);
     endpoint_vel_r = centraldiff(endpointright,120);
     endpoint_vel_magn_r = sqrt(sum(endpoint_vel_r.^2, 2));
     timev = 1:size(endpoint_vel,1);
     h39 = figure(39)
     clf
     hold on
     plot(timev, endpoint_vel_magn_l, 'b')
     hold on
     plot(timev, model_endpoint_vel_magn_l, 'r')
     plot(timev, sum(ang_contribs_l), 'm')
     plot(timev, endpoint_vel_magn_r, 'b:')
     plot(timev, model_endpoint_vel_magn_r, 'r:')
     plot(timev, sum(ang_contribs_r), 'm:')
     
     legend('Raw end point vel left', 'Model end point vel left', ...
	    'Sum of contrib left', ...
	    'Raw end point vel right', 'Model end point vel right', 'Sum of contrib right', ... 
	    'Location', 'Best')
     ylabel('Velocity [m/s]')
     xlabel('Frames')
    title(['Endpoint velocity, using ', epstr, ' as endpoint'],...
	   'Interpreter', 'None')
     
     
     
     % Find the impact event. 
     %events{1,1} = 'impact';
     %events{1,2} = find_impact(gm, states);
     
     
    % Path to write results to
    pth = fileparts(filestr);
    ndir = fullfile(pth,['results_',datestr(date,29)]);
    mkdir(ndir);
    
    [pth,mfname] = fileparts(filestr);

    fname_export = fullfile(ndir, [mfname,'_events.txt']);
    %export_values(ang_contribs, gm.gcnames(:,1), gm.gcnames(:,1), ...
%		  fname_export, events);
    
    % Write three tsv files. One with simulated markers and CoMs,
    % one with simulated markers and jcs, and one with real and
    % simulated markers.
    fname1 = fullfile(ndir, [mfname,'_endpoint.tsv']);
    fname2 = fullfile(ndir, [mfname,'_jc.tsv']);
    fname3 = fullfile(ndir, [mfname,'_residuals.tsv']);
    

    mnames = simnames;
    mnames = cat(1, mnames, comnames);
    fid = fopen(fname1, 'w');
    write3dtsv(putvalue(mdata{1}, 'MARKER_NAMES', mnames),...
	       cat(2, msim, com)*1000, fid);
    fclose(fid);

    mnames = simnames;
    mnames = cat(1, mnames, jcnames);
    fid = fopen(fname2, 'w');
    write3dtsv(putvalue(mdata{1}, 'MARKER_NAMES', mnames),...
	       cat(2, msim, jcs)*1000, fid);
    fclose(fid);

    mnames = getvalue(mdata{1}, 'MARKER_NAMES');
    mnames = cat(1, mnames, simnames);
    fid = fopen(fname3, 'w');
    write3dtsv(putvalue(mdata{1}, 'MARKER_NAMES', mnames),...
	       cat(2, mdata{2},msim)*1000, fid);
    fclose(fid);

    
    
  end  
end

  
  
    