%% Script interaction_moments_test
%
% Tests the interaction moments calculations, by simulating a
% simple test case: All angles constant zero, except sinsuoidal
% rotation of the trunk.

% Kjartan Halvorsen
% 2017-07-04

% The degrees of freedom to include in the analysis of interaction moments
% Based on the set of dofs that contributed almost all clubhead
% speed in the ISBS 2012 paper
%close all
%clear all

omega = 2; % rad/s of sinusoid
amplitude = pi/2; % Amplitude of trunk rotation
dt = 0.01; % Sampling time
N = 80;  % Number of samples
 
tt = (0:(N-1))*dt;
rotation = amplitude * sin(omega*tt);
angvel = amplitude*omega*cos(omega*tt); % Why cosine here?
angacc = -amplitude*omega^2*sin(omega*tt);

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

%datapth = ['C:\Users\fredrikt\Documents\MATLAB_OLD\EndpointContributions'];
datapth = ['/home/kjartan/Dropbox/SAS/Mobility'];

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
        startfr = getfield(getfield(startFrame,fps{fp,1}), trial);
        mdata{2} = mdata{2}(startfr:end,:);
        mdata{1,1}{1,2} = num2str(length(mdata{2}));       
    end
    
    % Create models.
    % The trial data are needed because the position of the club
    % with respect to either hand is taken from the first frame
    % (address) of the trial file.
    [gmleft, gmright, gmbase, gmclub] = ...
       build_im_models(refdata, mdata, bodymass);
  
    nstsleft = size(gmleft.gcnames, 1);
    
    rotationind_left = find(ismember(gmleft.gcnames(:,1), 'trunk rotation'));
    rotationind_base = find(ismember(gmbase.gcnames(:,1), 'trunk rotation'));
    
    flexind = find(ismember(gmleft.gcnames(:,1), 'left shoulder flexion'));
    shoulderaddind = find(ismember(gmleft.gcnames(:,1), 'left shoulder adduction'));
    elbowflexind = find(ismember(gmleft.gcnames(:,1), 'left elbow flexion'));

    % Apply sinusoidal position and velocity in trunk rotation.
    statesleft = zeros(nstsleft*2, N);
    statesleft(rotationind_left, :) = rotation;
    statesleft(nstsleft + rotationind_left, :) = angvel;
    
    % Constant 90 degrees in shoulder abduction
    statesleft(shoulderaddind, :) = -pi/2;
  
    % And constant 90 degrees in elbow flexion - should give more effect
    % form coriolis, less from acc
    %statesleft(elbowflexind, :) = pi/2;
    
    
    
    nstsright = size(gmright.gcnames, 1);
    statesright = zeros(nstsright*2, N);
    nstsbase = size(gmbase.gcnames, 1);
    statesbase = zeros(nstsbase*2, N);
    statesbase(rotationind_base, :) = rotation;
    statesbase(nstsbase + rotationind_base, :) = angvel;
    
    % Simulate models, generate trajectories of joint centers.
    % plot markers and check the residuals
    [nstsleft, nfrsleft] = size(statesleft);
    [msimleft, simnamesleft, objdleft, objnamesleft] = ...
        sim_model(gmleft, statesleft(1:nstsleft/2,:), 'jcs');
    [msimleft, simnamesleft, comleft, comnamesleft] = ...
        sim_model(gmleft, statesleft(1:nstsleft/2,:), 'CoM');
    [msimleft, simnamesleft, flexaxisleft, flexaxisnames] = ...
        sim_model(gmleft, statesleft(1:nstsleft/2,:), 'flexaxis');
     
    trunkind = find(ismember(objnamesleft, 'trunk_jc'));
    shoulderind = find(ismember(objnamesleft, 'left_upper_arm_jc'));
    elbowind = find(ismember(objnamesleft, 'left_forearm_jc'));
    trunkjoint = objdleft(:, ...
                             (trunkind-1)*3+1:(trunkind*3));
    shoulderjoint = objdleft(:, ...
                             (shoulderind-1)*3+1:(shoulderind*3));
    elbowjoint = objdleft(:, ...
                             (elbowind-1)*3+1:(elbowind*3));

    comUpperarmind = find(ismember(comnamesleft, 'left_upper_arm_CoM'));
    comForearmind = find(ismember(comnamesleft, 'left_forearm_CoM'));
    comHandind = find(ismember(comnamesleft, 'left_hand_CoM'));
    comUpperarmTimeseries = comleft(:, ...
                             (comUpperarmind-1)*3+1:(comUpperarmind*3));
    comForearmTimeseries = comleft(:, ...
                             (comForearmind-1)*3+1:(comForearmind*3));
    comHandTimeseries = comleft(:, ...
                             (comHandind-1)*3+1:(comHandind*3));
    
    elbowflexaxisind = find(ismember(flexaxisnames, 'left_forearm_flexaxis'));
    elbowflexaxis = flexaxisleft(:, ...
                                    (elbowflexaxisind-1)*3+1:(elbowflexaxisind*3));
    
                         
                         
                         
    rShoulderElbow = elbowjoint - shoulderjoint;
    rTrunkShoulder = shoulderjoint - trunkjoint;
    
 
    % The movement of the vector rTrunkShoulder in the direction of
    % the flexaxis should be constant
    flexaxis = gmleft.gcnames{rotationind_left, 3};
    tstMov = rTrunkShoulder * flexaxis;
    figure(1)
    clf
    plot(tstMov)
    title('This should be constant')
    
    
    % The angular velocity and -acceleration vectors for the left arm
    armVel = [flexaxis(1)*angvel', flexaxis(2)*angvel', flexaxis(3)*angvel']; % N x 3
    armAcc = [flexaxis(1)*angacc', flexaxis(2)*angacc', flexaxis(3)*angacc']; % N x 3
    
    % The acceleration of the elbow joint
    shoulderAcc = zeros(size(armAcc));
    elbowAcc = zeros(size(armAcc));
    elbowAccTangential = zeros(size(armAcc));
    for i=1:N
        shoulderAcc(i,:) = cross(armAcc(i,:), rTrunkShoulder(i,:)) ...
            + cross(armVel(i,:), cross(armVel(i,:), ...
                                                rTrunkShoulder(i,:)));
        elbowAccTangential(i,:) = cross(armAcc(i,:), rShoulderElbow(i,:));
        elbowAcc(i,:) = shoulderAcc(i,:) + cross(armAcc(i,:), rShoulderElbow(i,:)) ...
            + cross(armVel(i,:), cross(armVel(i,:), ...
                                                rShoulderElbow(i,:)));
    end

    % Now calculated by twice numerical differentiation.
    sfreq = 1.0/dt;
    elbowVelTest = centraldiff(elbowjoint, sfreq);
    elbowAccTest = centraldiff(elbowVelTest, sfreq);
    shoulderVelTest = centraldiff(shoulderjoint, sfreq);
    shoulderAccTest = centraldiff(shoulderVelTest, sfreq);
    
    figure(2)
    clf
    subplot(1,2,1)
    plot(elbowAcc, 'linewidth', 4)
    hold on
    plot(elbowAccTest)
    title(strcat('Elbow acc as  ', ...
                 '$ a \times r + \omega \times \omega \times r$', ...
                 ' (thick lines) and as $\frac{d^2}{dt^2} r$ ', ...
                 '(thin lines)'), 'interpreter', 'latex')
    
    subplot(1,2,2)
    plot(shoulderAcc, 'linewidth', 4)
    hold on
    plot(shoulderAccTest)
    title(strcat('Shoulder acc as  ', ...
                 '$ a \times r + \omega \times \omega \times r$', ...
                 ' (thick lines) and as $\frac{d^2}{dt^2} r$ ', ...
                 '(thin lines)'), 'interpreter', 'latex')
    
    
    
    % The moment of inertia of the forearm and hand wrt the elbow joint
    
    %                        pelv tru upper forearm
    IForearm = gmleft.inertia{ 2 }{2 }{ 2  }{1}(4:6, 4:6);
    IHand = gmleft.inertia{ 2 }{2 }{ 2  }{2}{1}(4:6, 4:6);
    mForearm = gmleft.mass{4};
    mHand = gmleft.mass{5};
    comForearm = gmleft.CoM{ 2 }{2 }{ 2  }{1}{2};
    comHand = gmleft.CoM{ 2 }{2 }{ 2  }{2}{1}{2};
    elbowRef = gmleft.jcs{2}{2}{2}{1}{2};

    [ILArm3, comLArm, mLArm] = combine_inertia(IForearm, comForearm, ...
                                              mForearm, ...
                                              IHand, comHand, ...
                                              mHand, elbowRef);
    
    
    % Moment of inertia wrt elbow flexion axis and rotation axis
    % Check against generalized inertia
    Mleft = generalized_manipulator_inertia(gmleft, statesleft);
    elbowflexaxislocal = gmleft.gcnames{elbowflexind, 3};
    ILArm = elbowflexaxislocal'*ILArm3*elbowflexaxislocal;
    Mleft(elbowflexind, elbowflexind, 1)
    
    figure(3)
    clf
    plot(reshape(Mleft(elbowflexind, elbowflexind, :), [N,1]))
    ylim([0.10, 0.11])
    title(sprintf('Should be constant and close to %f', ILArm))
    
    comLArmTimeseries = (1/(mHand+mForearm)) ... 
        * (mForearm*comForearmTimeseries + mHand* ...
           comHandTimeseries);
    
    
    %% Angle between elbowflexaxis and rotation axis
    angle = 180/pi * acos(elbowflexaxis*flexaxis);
    
    figure(11)
    clf
    plot(angle)
    
    %% The joint torque needed to produce the movement of the
    %% forarm and hand
    
    rLArm = norm(comLArm - elbowRef);
    rLArmTimeseries = comLArmTimeseries - elbowjoint;
    elbowAccMoment = zeros(N, 1);
    elbowAccMoment3 = zeros(N,3);
    elbowflexacc = zeros(N,1);
    
    for i=1:N
        elbowAccMoment3(i,:) = cross(rLArmTimeseries(i,:), ...
                                  -mLArm*elbowAcc(i,:));     % Since we are observing in accelerating frame of ref
        elbowAccMoment(i) = elbowAccMoment3(i,:)*elbowflexaxis(i,:)'; %% OBS: flexaxis not constant.
        elbowflexacc(i) = armAcc(i,:)*elbowflexaxis(i,:)'; % armAcc is the angular acc of the arm
    end
    
    %elbowAccMoment = norm(sqrt(sum(elbowAccTangential.^2, 2));
    %tauLArm = ILArm*angacc' + mLArm*rLArm)*elbowAccNormTangential.*sign(angacc');
    %tauLArm = ILArm*angacc' - mLArm*elbowAccMoment;
    tauLArmAngAcc = ILArm*elbowflexacc;
    tauLArmAccFrame = -elbowAccMoment;
    tauLArm = tauLArmAngAcc + tauLArmAccFrame; 
    % Negative signs since elbow flexaxis is opposite direction of
    % rotation axis
    


    %% The regular interaction moment calculations
    nfrs = size(statesleft, 2);
     nBaseStates = size(statesbase, 1)/2;
     nLeftArmStates = size(statesleft, 1)/2 - nBaseStates;
     nRightArmStates = size(statesright, 1)/2 - nBaseStates;
     nBS = nBaseStates;
     nLS = nLeftArmStates;
     nRS = nRightArmStates;
    
     Mbase = generalized_manipulator_inertia(gmbase, statesbase(1:nBS, :));
     Mleft = generalized_manipulator_inertia(gmleft, statesleft(1:nBS+nLS, :));
     Mright = generalized_manipulator_inertia(gmright, statesright(1:nBS+nRS, :));

     % Debug. OK!
     [twsBase, g0Base, MbBase, gcnamesBase] = flatten_km(gmbase);
     [twsLeft, g0Left, MbLeft, gcnamesLeft] = flatten_km(gmleft);
     [twsRight, g0Right, MbRight, gcnamesRight] = flatten_km(gmright);
     Mbase_flat = generalized_manipulator_inertia_flattened(MbBase, twsBase, g0Base, statesbase(1:nBS, 60));
     Mleft_flat = generalized_manipulator_inertia_flattened(MbLeft, twsLeft, g0Left, statesleft(1:nLS+nBS, 60));
     Mright_flat = generalized_manipulator_inertia_flattened(MbRight, twsRight, g0Right, statesright(1:nRS+nBS, 60));
     
     Mbase = Mbase(1:nBS, 1:nBS, :); % Kolla om inte har nBS rader från början
     Mleft = Mleft(1:nBS+nLS, 1:nBS+nLS, :);
     Mright = Mright(1:nBS+nRS, 1:nBS+nRS, :);
     
     
     % Combine
     ML_B = Mleft(1:nBS, 1:nBS,:);
     ML_L = Mleft(nBS+1:end,nBS+1:end,:);
     ML_BL = Mleft(1:nBS, nBS+1:end,:);
     MR_B = Mright(1:nBS, 1:nBS,:);
     MR_R = Mright(nBS+1:end,nBS+1:end,:);
     MR_BR = Mright(1:nBS, nBS+1:end,:);
 
     Mall = zeros(nBS+nLS+nRS, nBS+nLS+nRS, nfrs);
     Mall(1:nBS, 1:nBS, :) = Mbase+ML_B+MR_B;
    
     Mall(1:nBS, nBS+1:nBS+nLS, :) = ML_BL;
     Mall(nBS+1:nBS+nLS, 1:nBS,:) = permute(ML_BL, [2,1,3]);
     Mall(nBS+1:nBS+nLS, nBS+1:nBS+nLS,:) = ML_L;
 
     Mall(1:nBS, nBS+nLS+1:end, :) = MR_BR;
     Mall(nBS+1+nLS:end, 1:nBS,:) = permute(MR_BR, [2,1,3]);
     Mall(nBS+1+nLS:end, nBS+1+nLS:end,:) = MR_R;
     
     
     % And the joint accelerations
     freq=getvalue(mdata{1},'FREQUENCY');      % The sampling time
     if ischar(freq) freq=str2num(freq); end
     if (isempty(freq)) freq=240; end

     [nstsleft, nfrsleft] = size(statesleft);
     velLeft = statesleft(nstsleft/2+1:end,:);
     accLeft = centraldiff(velLeft', sfreq)'; % OBS sfreq is for this simulation!!!!
     
     [nstsright, nfrsright] = size(statesright);
     velRight = statesright(nstsright/2+1:end,:);
     accRight = centraldiff(velRight', sfreq)';
                       
     [nstsbase, nfrsbase] = size(statesbase);
     velBase = statesbase(nstsbase/2+1:end,:);
     accBase = centraldiff(velBase', sfreq)';

     accAll = cat(1, accBase, accLeft(nBS+1:end,:), accRight(nBS+1:end, :));
     
     % Calculate the interaction due to acceleration at other angles (degrees of
     % freedom). Remove the diagonwal, then multiply Mall with acc vector
     ntot = nBS+nLS+nRS;
     IM_acc = zeros(ntot, nfrs);
     for i=1:nfrs
         Mi = Mall(:,:,i); 
         Mi_nodiagonal = Mi - diag(diag(Mi));
         % IM_acc(:,i) = - Mi_nodiagonal*cat(1, accBase(:,i), accLeft(nBS+1:end,i), accRight(nBS+1:end,i));
         IM_acc(:,i) = - Mi_nodiagonal*accAll(:,i);
     end
     

     %% Calculate the interaction moments. OBS these are due to velocity only. Must be combined with IM_acc
     % im is the interaction moment for each degree of freedom, C
     % is the Coriolis matrix, i.e. 
     %  im = -C * \dot{q} 
     disp('Computing interaction moments for left arm')
     [IMleft, Cleft, dofnamesLeft] = interaction_moments(gmleft, statesleft); %, dofs2study);
     disp('Computing interaction moments for right arm')
     [IMright, Cright, dofnamesRight] = interaction_moments(gmright, statesright);%, dofs2study);
     disp('Computing interaction moments for hip and trunk')
     [IMbase, Cbase, dofnamesBase] = interaction_moments(gmbase, statesbase);%, dofs2study);

     %% Combine interaction terms 
     
     
     IMbaseLR = IMbase + IMleft(1:nBaseStates, :) ...
         + IMright(1:nBaseStates, :);

     IM_vel = cat(1, IMbaseLR, IMleft(nBaseStates+1:end,:), ...
         IMright(nBaseStates+1:end, :));
     IMall =  IM_vel + IM_acc;
     
 
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
             

     %% Now the induced accelerations. Since we have taken care of the off-diagonal inertial terms these should not be included
     indAcc = zeros(ntot, nfrs);
     for i = 1:nfrs
         indAcc(:,i) = IMall(:, i)./diag(Mall(:,:,i));
     end
     
     
     
     %% Verify the calculated interaction moment of the elbow
     %% flexion
     tauLArm2 = IMall(elbowflexind, :)';
     tauLArm2Vel = IM_vel(elbowflexind, :)';
     tauLArm2Acc = IM_acc(elbowflexind, :)';
     
     % Plot interaction moment in elbow flexion, and the computed joint torque
     figure(4)
     clf
     plot(tauLArm2, 'linewidth', 3)
     hold on
     plot(tauLArm, 'linewidth', 2)
     plot(tauLArm2Vel)
     plot(tauLArm2Acc)
     plot(tauLArmAccFrame)
     plot(tauLArmAngAcc)
     legend('Interaction moment', 'elbow joint torque', ...
         'Interaction moment, velocity dependent', ...
         'Interaction moment, acc dependent', ...
         'joint torque due to accelerating frame fixed in elbow', ...
         'joint torque due to angular acceleration')
     
    

     
     
     
     if debug
         dofnames = [dofnamesLeft; dofnamesRight(nBaseStates+1:end)];
         % Plot the interaction moments
         figure(5)
         clf
         for dof=1:nLeftArmStates
            subplot(ceil(nLeftArmStates/2),2, dof)
            plot(IMleft(nBaseStates+dof,:)')
            title(dofnamesLeft{nBaseStates+dof})
         end
         
%          subplot(121)
%          plot(IMleft(nBaseStates+1:end,:)')
%          title('Interaction moments, left arm')
%          legend(dofnamesLeft{nBaseStates+1:end})
%          subplot(122)
%          plot(IMright(nBaseStates+1:end,:)')
%          legend(dofnamesRight{nBaseStates+1:end})
%          title('Interaction moments, right arm')
%          
         figure(6)
         clf
         for dof=1:nRightArmStates
            subplot(ceil(nRightArmStates/2),2, dof)
            plot(IMright(nBaseStates+dof,:)')
            title(dofnamesRight{nBaseStates+dof})
         end
%         plot(IMbaseLR')
%         title('Interaction moments, hip and trunk')
%         legend(dofnamesBase)
      
         figure(7)
         clf
         for dof=1:nLeftArmStates
            subplot(ceil(nLeftArmStates/2),2, dof)
            plot(accLeft(nBaseStates+dof,:)')
            hold on
            plot(indAcc(nBaseStates+dof,:)')
            title(dofnamesLeft{nBaseStates+dof})
         end
         
         figure(8)
         clf
         for dof=1:nRightArmStates
            subplot(ceil(nRightArmStates/2),2, dof)
            plot(accRight(nBaseStates+dof,:)')
            hold on
            plot(indAcc(nBaseStates+nLeftArmStates+dof,:)')
            title(dofnamesRight{nBaseStates+dof})
         end
         
        
     end

    
  end
end

