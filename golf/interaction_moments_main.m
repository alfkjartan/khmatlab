% Script interaction_moments_main
% Tracks joint angles and computes interaction moments 

% Based on track_main_FT

% Kjartan Halvorsen
% 2017-04-19

% The degrees of freedom to include in the analysis of interaction moments
% Based on the set of dofs that contributed almost all clubhead
% speed in the ISBS 2012 paper
dofs2study = {'pelvis obliquety', ...
              'pelvis rotation', ...
              'trunk rotation', ...
              'left shoulder flexion', ...
              'left shoulder abduction', ...
              'left shoulder rotation', ...
              'left elbow rotation', ...
              'left wrist flexion', ...
              'left wrist abduction', ...
              'right shoulder flexion', ...
              'right shoulder abduction', ...
              'right shoulder rotation', ...
              'right elbow flexion', ...
              'right wrist flexion', ...
              'right wrist abduction'};
              

debug = 1;

filterbandwidth = 1250;  % Scales the process noise covariance
                      % matrix. Increase the value for higher
                      % bandwidth (less smoothing effect), lower
                      % for more smoothing.

datapth = ['C:\Users\fredrikt\Documents\MATLAB\EndpointContributions'];
datapth = ['/home/kjartan/Dropbox/SAS/Mobility'];

% Initials, folder name and reference mat-file
% fps = {'AH', 'AH', 'wedgemodel_ah.mat'
%        'AW', 'AW', 'wedgemodel_aw.mat'
%        'DK', 'DK', 'wedgemodel_dk.mat'
%        'EE', 'EE', 'wedgemodel_ee.mat'
%        'GN', 'GN', 'wedgemodel_gn.mat'
%        'JJ', 'JJ', 'wedgemodel_jj.mat'
%        'KH', 'KH', 'wedgemodel_kh.mat'
%        'MBE', 'MBE', 'wedgemodel_mbe.mat'
%        'SN', 'SN', 'wedgemodel_sn.mat'
%        'SP', 'SP', 'wedgemodel_sp.mat'};
   
% Initials, folder name, reference mat-file and body mass
fps = {'AH', 'AH', 'wedgemodel_AH.mat', 80
       'AW', 'AW', 'wedgemodel_aw.mat', 80
       'DK', 'DK', 'wedgemodel_dk.mat', 80
       'EE', 'EE', 'wedgemodel_ee.mat', 80
       'GN', 'GN', 'wedgemodel_gn.mat', 80
       'JJ', 'JJ', 'wedgemodel_jj.mat', 80
       'KH', 'KH', 'wedgemodel_kh.mat', 80
       'MBE', 'MBE', 'wedgemodel_mbe.mat', 80
       'SN', 'SN', 'wedgemodel_sn.mat', 80
       'SP', 'SP', 'wedgemodel_sp.mat', 80
       'DP', 'DP', 'wedgemodel_dp.mat', 80
       'FP', 'FP', 'wedgemodel_fp.mat', 80
       'LJ', 'LJ', 'wedgemodel_lj.mat', 80
       'MW', 'MW', 'wedgemodel_mw.mat', 80
       'AL', 'AL', 'wedgemodel_al.mat', 80
       'EA', 'EA', 'wedgemodel_ea.mat', 80
       'FB', 'FB', 'wedgemodel_fb.mat', 80
       'MB', 'MB', 'wedgemodel_mb.mat', 80
       'NB', 'NB', 'wedgemodel_nb.mat', 80
       'SS', 'SS', 'wedgemodel_ss.mat', 80};

trials = {'File251.c3d'
	  'File252.c3d'
	  'File253.c3d'
	  'File551.c3d'
	  'File552.c3d'
	  'File553.c3d'
	  'Filefull1.c3d'
	  'Filefull2.c3d'
	  'Filefull3.c3d'};
  
%   trials = {'File251.c3d'
% 	  'File252.c3d'
% 	  'File253.c3d'
% 	  'File401.c3d'
% 	  'File402.c3d'
% 	  'File403.c3d'
% 	  'File551.c3d'
% 	  'File552.c3d'
% 	  'File553.c3d'
% 	  'File701.c3d'
% 	  'File702.c3d'
% 	  'File703.c3d'
% 	  'Filefull1.c3d'
% 	  'Filefull2.c3d'
% 	  'Filefull3.c3d'};

useCustomStartFrame = 1;

startFrame.DP.File251=10;
startFrame.DP.File252=40;
startFrame.DP.File253=200;
startFrame.DP.File401=150;
startFrame.DP.File402=180;
startFrame.DP.File403=230;
startFrame.DP.File551=290;
startFrame.DP.File552=1;
startFrame.DP.File553=20;
startFrame.DP.File701=100;
startFrame.DP.File702=290;
startFrame.DP.File703=190;
startFrame.DP.Filefull1=290;
startFrame.DP.Filefull2=400;
startFrame.DP.Filefull3=410;

startFrame.FP.File251=90;
startFrame.FP.File252=350;
startFrame.FP.File253=320;
startFrame.FP.File401=1;
startFrame.FP.File402=150;
startFrame.FP.File403=210;
startFrame.FP.File551=70;
startFrame.FP.File552=150;
startFrame.FP.File553=100;
startFrame.FP.File701=240;
startFrame.FP.File702=140;
startFrame.FP.File703=220;
startFrame.FP.Filefull1=330;
startFrame.FP.Filefull2=250;
startFrame.FP.Filefull3=250;

startFrame.LJ.File251=210;
startFrame.LJ.File252=210;
startFrame.LJ.File253=80;
startFrame.LJ.File401=190;
startFrame.LJ.File402=250;
startFrame.LJ.File403=50;
startFrame.LJ.File551=200;
startFrame.LJ.File552=10;
startFrame.LJ.File553=50;
startFrame.LJ.File701=150;
startFrame.LJ.File702=180;
startFrame.LJ.File703=250;
startFrame.LJ.Filefull1=100;
startFrame.LJ.Filefull2=200;
startFrame.LJ.Filefull3=200;

startFrame.MW.File251=470;
startFrame.MW.File252=370;
startFrame.MW.File253=320;
startFrame.MW.File401=370;
startFrame.MW.File402=940;
startFrame.MW.File403=290;
startFrame.MW.File551=500;
startFrame.MW.File552=660;
startFrame.MW.File553=340;
startFrame.MW.File701=540;
startFrame.MW.File702=500;
startFrame.MW.File703=330;
startFrame.MW.Filefull1=450;
startFrame.MW.Filefull2=400;
startFrame.MW.Filefull3=350;

startFrame.MH.File251=40;
startFrame.MH.File401=120;
startFrame.MH.File551=190;
startFrame.MH.File701=350;
startFrame.MH.Filefull1=270;

startFrame.SS.File251=1;
startFrame.SS.File252=1;
startFrame.SS.File401=20;
startFrame.SS.File402=50;
startFrame.SS.File551=30;
startFrame.SS.File552=1350;
startFrame.SS.File701=1440;
startFrame.SS.File702=60;
startFrame.SS.Filefull1=20;
startFrame.SS.Filefull2=100;

startFrame.AL.File251=265;
startFrame.AL.File252=1060;
startFrame.AL.File253=570;
startFrame.AL.File401=300;
startFrame.AL.File402=1;
startFrame.AL.File403=100;
startFrame.AL.File551=40;
startFrame.AL.File552=215;
startFrame.AL.File553=160;
startFrame.AL.File701=125;
startFrame.AL.File702=60;
startFrame.AL.File703=1;
startFrame.AL.Filefull1=435;
startFrame.AL.Filefull2=840;
startFrame.AL.Filefull3=160;

startFrame.EA.File251=260;
startFrame.EA.File252=520;
startFrame.EA.File253=340;
startFrame.EA.File401=300;
startFrame.EA.File402=270;
startFrame.EA.File403=270;
startFrame.EA.File551=240;
startFrame.EA.File552=350;
startFrame.EA.File553=330;
startFrame.EA.File701=150;
startFrame.EA.File702=260;
startFrame.EA.File703=130;
startFrame.EA.Filefull1=435;
startFrame.EA.Filefull2=320;
startFrame.EA.Filefull3=310;

startFrame.FB.File251=470;
startFrame.FB.File252=270;
startFrame.FB.File253=390;
startFrame.FB.File401=120;
startFrame.FB.File402=290;
startFrame.FB.File403=210;
startFrame.FB.File551=120;
startFrame.FB.File552=290;
startFrame.FB.File553=300;
startFrame.FB.File701=130;
startFrame.FB.File702=100;
startFrame.FB.File703=280;
startFrame.FB.Filefull1=290;
startFrame.FB.Filefull2=350;
startFrame.FB.Filefull3=300;

startFrame.MB.File251=420;
startFrame.MB.File252=380;
startFrame.MB.File253=470;
startFrame.MB.File401=650;
startFrame.MB.File402=70;
startFrame.MB.File403=230;
startFrame.MB.File551=30;
startFrame.MB.File552=520;
startFrame.MB.File553=750;
startFrame.MB.File701=370;
startFrame.MB.File702=310;
startFrame.MB.File703=280;
startFrame.MB.Filefull1=360;
startFrame.MB.Filefull2=400;
startFrame.MB.Filefull3=390;

startFrame.NB.File251=250;
startFrame.NB.File252=200;
startFrame.NB.File253=250;
startFrame.NB.File401=320;
startFrame.NB.File402=320;
startFrame.NB.File403=210;
startFrame.NB.File551=390;
startFrame.NB.File552=60;
startFrame.NB.File553=200;
startFrame.NB.File701=250;
startFrame.NB.File702=250;
startFrame.NB.File703=170;
startFrame.NB.Filefull1=320;
startFrame.NB.Filefull2=230;
startFrame.NB.Filefull3=260;

startFrame.AH.File251=240;
startFrame.AH.File252=400;
startFrame.AH.File253=315;
startFrame.AH.File551=450;
startFrame.AH.File552=410;
startFrame.AH.File553=290;
startFrame.AH.Filefull1=430;
startFrame.AH.Filefull2=405;
startFrame.AH.Filefull3=860;

startFrame.AW.File251=640;
startFrame.AW.File252=120;
startFrame.AW.File253=1230;
startFrame.AW.File551=270;
startFrame.AW.File552=590;
startFrame.AW.File553=620;
startFrame.AW.Filefull1=810;
startFrame.AW.Filefull2=820;
startFrame.AW.Filefull3=995;

startFrame.DK.File251=45;
startFrame.DK.File252=55;
startFrame.DK.File253=365;
startFrame.DK.File551=50;
startFrame.DK.File552=125;
startFrame.DK.File553=120;
startFrame.DK.Filefull1=110;
startFrame.DK.Filefull2=270;
startFrame.DK.Filefull3=160;

startFrame.EE.File251=610;
startFrame.EE.File252=790;
startFrame.EE.File253=605;
startFrame.EE.File551=550;
startFrame.EE.File552=680;
startFrame.EE.File553=840;
startFrame.EE.Filefull1=300;
startFrame.EE.Filefull2=670;
startFrame.EE.Filefull3=650;

startFrame.GN.File251=155;
startFrame.GN.File252=270;
startFrame.GN.File253=75;
startFrame.GN.File551=280;
startFrame.GN.File552=240;
startFrame.GN.File553=1110;
startFrame.GN.Filefull1=440;
startFrame.GN.Filefull2=315;
% startFrame.GN.Filefull3=315;
startFrame.GN.Filefull3=510;

startFrame.JJ.File251=730;
startFrame.JJ.File252=660;
startFrame.JJ.File253=350;
startFrame.JJ.File551=730;
startFrame.JJ.File552=580;
startFrame.JJ.File553=1200;
startFrame.JJ.Filefull1=900;
startFrame.JJ.Filefull2=310;
startFrame.JJ.Filefull3=570;

startFrame.KH.File251=10;
startFrame.KH.File252=30;
startFrame.KH.File253=10;
startFrame.KH.File551=40;
startFrame.KH.File552=90;
startFrame.KH.File553=110;
startFrame.KH.Filefull1=145;
startFrame.KH.Filefull2=40;
startFrame.KH.Filefull3=320;

startFrame.MBE.File251=260;
startFrame.MBE.File252=20;
startFrame.MBE.File253=110;
startFrame.MBE.File551=550;
startFrame.MBE.File552=105;
startFrame.MBE.File553=145;
startFrame.MBE.Filefull1=295;
startFrame.MBE.Filefull2=125;
startFrame.MBE.Filefull3=410;

startFrame.SN.File251=620;
startFrame.SN.File252=770;
startFrame.SN.File253=620;
startFrame.SN.File551=540;
startFrame.SN.File552=690;
startFrame.SN.File553=830;
startFrame.SN.Filefull1=310;
startFrame.SN.Filefull2=680;
% startFrame.SN.Filefull3=680;
startFrame.SN.Filefull3=660;

startFrame.SP.File251=450;
startFrame.SP.File252=1;
startFrame.SP.File253=170;
startFrame.SP.File551=160;
startFrame.SP.File552=660;
startFrame.SP.File553=60;
% startFrame.SP.Filefull1=410;
startFrame.SP.Filefull1=1;
startFrame.SP.Filefull2=410;
% startFrame.SP.Filefull3=410;
startFrame.SP.Filefull3=100;

usefps = (2); % Select all or part of data to process
%usetrials = (13:14);
usetrials = (7);

% Close all open plot windows
close all

% Names of markers to plot must match exactly with the names in the
% c3d file.
 markers2plot = {'L_HAND_1'
                 'L_HAND_2'
                 'L_HAND_3'
                 'Club_1'
                 'Club_2'
                 'Club_3'
                 'ClubCoM'};
 for i=1:length(markers2plot)
   markers2plot{i,2} = figure('Name', markers2plot{i,1});
 end


% Names of joint angles to plot must match exactly with the names in the
% build_golf_model.m file.
angles2plot = {'pelvis x'
             'pelvis y'
             'pelvis z'
             'pelvis tilt'
             'pelvis obliqueity'
             'pelvis rotation'
             'trunk tilt'
             'trunk obliquety'
             'trunk rotation'
             'left shoulder x'
             'left shoulder y'
             'left shoulder z'
             'left shoulder flexion'
             'left shoulder abduction'
             'left shoulder rotation'
             'left elbow flexion'
             'left elbow rotation'
             'left wrist flexion'
             'left wrist abduction'
             'club_l x'
             'club_l y'
             'club_l z'
             'club_l tilt'
             'club_l yaw'
             'club_l rotation'
             'right shoulder x'
             'right shoulder y'
             'right shoulder z'
             'right shoulder flexion'
             'right shoulder abduction'
             'right shoulder rotation'
             'right elbow flexion'
             'right elbow rotation'
             'right wrist flexion'
             'right wrist abduction'
             'club_r x'
             'club_r y'
             'club_r z'
             'club_r tilt'
             'club_r yaw'
             'club_r rotation'};

angles2plot = dofs2study;

% for i=1:length(angles2plot)
%   angles2plot{i,2} = figure('Name', angles2plot{i,1});
% end

% Names of joint angles to plot must match exactly with the names in the
% build_golf_model.m file.
% anglecontribs2plot = {'pelvis rotation'
% 		    'trunk rotation'
% 		    'left shoulder flexion'
%             'left shoulder abduction'
%             'left shoulder rotation'
% 		    'left elbow flexion'
%             'left elbow rotation'
% 		    'left wrist flexion'
% 		    'left wrist abduction'};
% 
% for i=1:length(anglecontribs2plot)
%   anglecontribs2plot{i,2} = figure('Name', ...
% 				   [anglecontribs2plot{i,1},' contrib']);
% end

usemdata = [76:84, 154:162, 61:69, 52:60, 10:18, 130:138, 121:129, 10:18
            76:84, 157:165, 61:69, 52:60, 10:18, 130:138, 121:129, 10:18
            76:84, 160:168, 61:69, 52:60, 10:18, 130:138, 121:129, 10:18
            88:96, 166:174, 67:75, 58:66, 10:18, 142:150, 133:141, 10:18
            91:99, 172:180, 70:78, 61:69, 10:18, 148:156, 139:147, 10:18
            100:108, 181:189, 67:75, 58:66, 10:18, 154:162, 145:153, 10:18];
        
RMSE = zeros(length(usefps),(length(usemdata)/3));

% Now process the data
for fp=usefps
  
  bodymass = fps{fp,4};
  refdata = load(fullfile(datapth,fps{fp,2},fps{fp,3}));
    
  
  for tr=usetrials
    
    % Read the marker data
    filestr = fullfile(datapth,fps{fp,2},trials{tr});
    mdata = openmocapfile('', filestr);

    if useCustomStartFrame
        trial = trials{tr}(1:end-4);
        startfr = getfield(getfield(startFrame,fps{fp,1}), trial)
        mdata{2} = mdata{2}(startfr:end,:);
        mdata{1,1}{1,2} = num2str(length(mdata{2}));       
    end
    
    % Create models.
    % The trial data are needed because the position of the club
    % with respect to either hand is taken from the first frame
    % (address) of the trial file.
    [gmleft, gmright, gmbase, gmclub] = ...
       build_im_models(refdata, mdata, bodymass);

    
    [statesleft, dataframesleft] = track_golf_model(gmleft, mdata, filterbandwidth);
    [statesright, dataframesright] = track_golf_model(gmright, mdata, filterbandwidth);
    [statesbase, dataframesbase] = track_golf_model(gmbase, mdata, filterbandwidth);
    
     if debug
         % Simulate models, generate trajectories of joint centers.
         % plot markers and check the residuals
         [nstsleft, nfrsleft] = size(statesleft);
         [msimleft, simnamesleft, objdleft, objnamesleft] = ...
             sim_model(gmleft, statesleft(1:nstsleft/2,:), 'objectcenter');

         [nstsright, nfrsright] = size(statesright);
         [msimright, simnamesright, objdright, objnamesright] = ...
             sim_model(gmright, statesright(1:nstsright/2,:), 'objectcenter');

         % Plot results

         y_observations=extractmarkers(mdata, simnamesleft);
         y_observations(find(y_observations==0)) = NaN;
         mse = plotmarkers(y_observations(find(dataframesleft),:), ...
                           simnamesleft,...
                           msimleft, simnamesleft, markers2plot, {'data', 'model'});

         %keyboard
         y_observations=extractmarkers(mdata, simnamesright);
         y_observations(find(y_observations==0)) = NaN;
         mse = plotmarkers(y_observations(find(dataframesright),:), ...
                           simnamesright,...
                           msimright, simnamesright, markers2plot, {'data', 'model'});
         
         convert_radians = 1;
         %plotangles(statesleft, gmleft.gcnames(:,1), angles2plot, convert_radians);

     end


     % Calculate the interaction moments. 
     % im is the interaction moment for each degree of freedom, C
     % is the Coriolis matrix, i.e. 
     %  im = C * \dot{q} 
     [IMleft, Cleft, dofnamesLeft] = interaction_moments(gmleft, statesleft, dofs2study);
     [IMright, Cright, dofnamesRight] = interaction_moments(gmright, statesright, dofs2study);
     [IMbase, Cbase, dofnamesBase] = interaction_moments(gmbase, statesbase, dofs2study);

     % Combine interaction terms for the generalized coordinates of
     % the base model
     nfrs = size(IMbase, 2);
     nBaseStates = size(IMbase, 1);
     nLeftArmStates = size(IMleft, 1) - nBaseStates;
     nRightArmStates = size(IMright, 1) - nBaseStates;
     
     IMbaseLR = IMbase + IMleft(1:nBaseStates, :) ...
         + IMright(1:nBaseStates, :);

     % When calculating the Coriolis matrix, we assume the vector
     % of joint velocities to be [\dot{q}_B, \dot{q}_L, \dot{q}_R]
     CbaseLR = zeros(nBaseStates, ...
                     nBaseStates+nLeftArmStates+nRightArmStates, ... 
                     nfrs); 
     CbaseLR(:,1:nBaseStates, :) = Cbase ...
         + Cleft(1:nBaseStates, 1:nBaseStates,:) ...
         + Cright(1:nBaseStates, 1:nBaseStates,:);
     CbaseLR(:, nBaseStates+1:nBaseStates+nLeftArmStates,:) = ...
         Cleft(1:nBaseStates, nBaseStates+1:end, :);
     CbaseLR(:, nBaseStates+nLeftArmStates+1:end, :) = ...
             Cright(1:nBaseStates, nBaseStates+1:end, :);
             

     if debug
         % Plot the interaction moments
         figure(1)
         clf
         subplot(121)
         plot(IMleft(nBaseStates+1:end,:)')
         title('Interaction moments, left arm')
         legend(dofnamesLeft{nBaseStates+1:end})
         subplot(122)
         plot(IMright(nBaseStates+1:end,:)')
         legend(dofnamesRight{nBaseStates+1:end})
         title('Interaction moments, right arm')
         
         figure(2)
         clf
         plot(IMbaseLR')
         title('Interaction moments, hip and trunk')
         legend(dofnamesBase)
         
     end
  end
end

  
  
    