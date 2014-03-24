function [tm, endpointstr] = build_golf_model_w_inertia(refdata,trialdata, bodymass, clubshaftmass, clubheadmass)
%  gm = build_golf_model_w_inertia(refdata, trialdata, bodymass, clubshaftmass, clubheadmass)
%
% Returns a kinematic model for the golf study.
% The model contains 9 segments:
%    pelvis, trunk, upperarms, underarms, hands, club
% The pelvis is the root segment.
%
% Input
%    refata    ->  Reference marker data. Must be struct exported
%                  from visual3D
%    trialdata ->  marker data from a swing trial. The start of the
%                  file should correspond to address position. Most
%                  importantly, that the grip is the same as the
%                  one used througout the swing.Must be struct exported
%                  from visual3D
%    bodymass  ->  Total mass of the subject
%    clubshaftmass ->  Mass of the shaft of the club
%    clubheadmass ->  Mass of the shaft of the club
%                  To make a reasonable approximation, the club is divided into bodies: 
%                  shaft and head.
%                  Data on the respective weight and shape of the two parts is approximated:
%                  shaft: assumed a tapered tube with 10mm radius at top and
%                  3mm at the hose.
%                  head: Homogenous ellipsoid.
% Output
%    tm        <-  struct. Contains the fields 
%        tws       <-  nested cell array of twists
%        gcnames   <-  array with names for the generalized
%                      coordinates.
%        p0        <-  nested cell with reference marker positions.
%        jcs       <-  nested cell with joint centra 
%        CoM       <-  nested cell with center of mass 
%    endpointstr   <-  If not empty, contains the name of the
%                      marker used as endpoint
%

% Kjartan Halvorsen
% 2009-06-26
% Based on build_tennis_model.m from  2003-09-16
  
% Revisions
% 2010-03-21   Added output endpointstr 
% 2010-06-20   Added 3dof for the right arm (FT)
% 2010-06-20   Added relative peak acceleration for each degree of freedom
%              during the downswing phase (FT)
% 2013-06-06   Added inertial parameters to the model. These are based on de Leva, 
%              J Biomech 1996 (KH)


% Find suitable frame to use as address position

club_1_tr = extractmarkers(trialdata, 'ClubCoM');
[impact, impact_fit, pquad, dist2address, max_before_backsw] ...
= find_impact_from_point(club_1_tr);

trialdata{2} = trialdata{2}(max_before_backsw:max_before_backsw+1,:);

c7 = mmyextractmeanmarkers(refdata, 'C7');
ij = mmyextractmeanmarkers(refdata, 'IJ');
shoulder_l = mmyextractmeanmarkers(refdata, 'L_Acromion');
shoulder_r = mmyextractmeanmarkers(refdata, 'R_Acromion');
ghjc_l = mmyextractmeanmarkers(refdata, 'Wrt_LShoulder');
ghjc_r = mmyextractmeanmarkers(refdata, 'Wrt_RShoulder');
elbow_lat_l = mmyextractmeanmarkers(refdata, 'L_Elbow_lateral');
elbow_med_l = mmyextractmeanmarkers(refdata, 'L_Elbow_medial');
elbow_lat_r = mmyextractmeanmarkers(refdata, 'R_Elbow_lateral');
elbow_med_r = mmyextractmeanmarkers(refdata, 'R_Elbow_medial');
wrist_radial_l = mmyextractmeanmarkers(refdata, 'L_Radial_wrist');
wrist_ulnar_l = mmyextractmeanmarkers(refdata, 'L_Ulnar_wrist');
wrist_radial_r = mmyextractmeanmarkers(refdata, 'R_Radial_wrist');
wrist_ulnar_r = mmyextractmeanmarkers(refdata, 'R_Ulnar_wrist');
asis_l = mmyextractmeanmarkers(refdata, 'L_ASIS');
asis_r = mmyextractmeanmarkers(refdata, 'R_ASIS');
psis_l = mmyextractmeanmarkers(refdata, 'L_PSIS');
psis_r = mmyextractmeanmarkers(refdata, 'R_PSIS');
t8 = mmyextractmeanmarkers(refdata, 'T8');
pelvis_1 = mmyextractmeanmarkers(refdata, 'Pelvis_1');
pelvis_2 = mmyextractmeanmarkers(refdata, 'Pelvis_2');
pelvis_3 = mmyextractmeanmarkers(refdata, 'Pelvis_3');
ut_1 = mmyextractmeanmarkers(refdata, 'Upper_Torso_1');
ut_2 = mmyextractmeanmarkers(refdata, 'Upper_Torso_2');
ut_3 = mmyextractmeanmarkers(refdata, 'Upper_Torso_3');
uarm_1_l = mmyextractmeanmarkers(refdata, 'L_Upper_Arm_1');
uarm_2_l = mmyextractmeanmarkers(refdata, 'L_Upper_Arm_2');
uarm_3_l = mmyextractmeanmarkers(refdata, 'L_Upper_Arm_3');
uarm_1_r = mmyextractmeanmarkers(refdata, 'R_Upper_Arm_1');
uarm_2_r = mmyextractmeanmarkers(refdata, 'R_Upper_Arm_2');
uarm_3_r = mmyextractmeanmarkers(refdata, 'R_Upper_Arm_3');
hand_1_l = mmyextractmeanmarkers(refdata, 'L_Hand_1');
hand_2_l = mmyextractmeanmarkers(refdata, 'L_Hand_2');
hand_3_l = mmyextractmeanmarkers(refdata, 'L_Hand_3');
hand_1_r = mmyextractmeanmarkers(refdata, 'R_Hand_1');
hand_2_r = mmyextractmeanmarkers(refdata, 'R_Hand_2');
hand_3_r = mmyextractmeanmarkers(refdata, 'R_Hand_3');
hand_1_l_trial = mmyextractmeanmarkers(trialdata, 'L_Hand_1');
hand_2_l_trial = mmyextractmeanmarkers(trialdata, 'L_Hand_2');
hand_3_l_trial = mmyextractmeanmarkers(trialdata, 'L_Hand_3');
hand_1_r_trial = mmyextractmeanmarkers(trialdata, 'R_Hand_1');
hand_2_r_trial = mmyextractmeanmarkers(trialdata, 'R_Hand_2');
hand_3_r_trial = mmyextractmeanmarkers(trialdata, 'R_Hand_3');
mp_2_l = mmyextractmeanmarkers(refdata, 'L_2nd_MP_joint');
mp_5_l = mmyextractmeanmarkers(refdata, 'L_5th_MP_joint');
mp_2_r = mmyextractmeanmarkers(refdata, 'R_2nd_MP_joint');
mp_5_r = mmyextractmeanmarkers(refdata, 'R_5th_MP_joint');
grip_top = mmyextractmeanmarkers(refdata, 'TopOfHandle');
heel_bottom_grove = mmyextractmeanmarkers(refdata, 'BottomGroveHeel');
toe_bottom_grove = mmyextractmeanmarkers(refdata, 'BottomGroveToe');
toe_top_grove = mmyextractmeanmarkers(refdata, 'TopGroveToe');
club_1 = mmyextractmeanmarkers(refdata, 'Club_1');
club_2 = mmyextractmeanmarkers(refdata, 'Club_2');
club_3 = mmyextractmeanmarkers(refdata, 'Club_3');
club_1_trial = mmyextractmeanmarkers(trialdata, 'Club_1');
club_2_trial = mmyextractmeanmarkers(trialdata, 'Club_2');
club_3_trial = mmyextractmeanmarkers(trialdata, 'Club_3');
%midhands_1 = mmyextractmeanmarkers(refdata, 'MidHands_1');
%midhands_2 = mmyextractmeanmarkers(refdata, 'MidHands_2');



% Transform the club markers.
g_cl_trial_ref = ...
    soder(cat(1, ...
	      cat(2, club_1', club_2', club_3'),...
	      cat(2, club_1_trial', club_2_trial', ...
		  club_3_trial')));
g_rh_ref_trial = ...
    soder(cat(1, ...
	      cat(2, hand_1_r_trial', hand_2_r_trial', ...
		  hand_3_r_trial'),...
	      cat(2, hand_1_r', hand_2_r', hand_3_r')));
g_lh_ref_trial = ...
    soder(cat(1, ...
	      cat(2, hand_1_l_trial', hand_2_l_trial', ...
		  hand_3_l_trial'),...
	      cat(2, hand_1_l', hand_2_l', hand_3_l')));

club_1_r = g_rh_ref_trial*g_cl_trial_ref * cat(1, club_1, 1);
club_1_r(4)=[];
club_2_r = g_rh_ref_trial*g_cl_trial_ref * cat(1, club_2, 1);
club_2_r(4)=[];
club_3_r = g_rh_ref_trial*g_cl_trial_ref * cat(1, club_3, 1);
club_3_r(4)=[];
grip_top_r = g_rh_ref_trial*g_cl_trial_ref * cat(1, grip_top, 1);
grip_top_r(4)=[];
heel_bottom_grove_r = g_rh_ref_trial*g_cl_trial_ref * cat(1, heel_bottom_grove, 1);
heel_bottom_grove_r(4)=[];
toe_bottom_grove_r = g_rh_ref_trial*g_cl_trial_ref * cat(1, toe_bottom_grove, 1);
toe_bottom_grove_r(4)=[];
toe_top_grove_r = g_rh_ref_trial*g_cl_trial_ref * cat(1, toe_top_grove, 1);
toe_top_grove_r(4)=[];

club_1_l = g_lh_ref_trial*g_cl_trial_ref * cat(1, club_1, 1);
club_1_l(4)=[];
club_2_l = g_lh_ref_trial*g_cl_trial_ref * cat(1, club_2, 1);
club_2_l(4)=[];
club_3_l = g_lh_ref_trial*g_cl_trial_ref * cat(1, club_3, 1);
club_3_l(4)=[];
grip_top_l = g_lh_ref_trial*g_cl_trial_ref * cat(1, grip_top, 1);
grip_top_l(4)=[];
heel_bottom_grove_l = g_lh_ref_trial*g_cl_trial_ref * cat(1, heel_bottom_grove, 1);
heel_bottom_grove_l(4)=[];
toe_bottom_grove_l = g_lh_ref_trial*g_cl_trial_ref * cat(1, toe_bottom_grove, 1);
toe_bottom_grove_l(4)=[];
toe_top_grove_l = g_lh_ref_trial*g_cl_trial_ref * cat(1, toe_top_grove, 1);
toe_top_grove_l(4)=[];


		       
% Landmarks and anatomical directions
midtrunk = mean(cat(2, shoulder_l,shoulder_r, ...
		    c7, ij), 2);
midpelvis = mean(cat(2, asis_l, asis_r, psis_r, psis_l), 2); 

e_IS = midtrunk - midpelvis; % Inferior-superior direction 
e_IS = e_IS / norm(e_IS);
e_LR = shoulder_r - shoulder_l; % Left-right, local X
e_LR = e_LR - (e_LR'*e_IS)*e_IS;
e_LR = e_LR/norm(e_LR);
e_PA = cross(e_IS, e_LR); % Posterior-anterior, local Y.


% Elbow joint centers
ejc_l = 0.5*elbow_lat_l + 0.5*elbow_med_l;
ejc_r = 0.5*elbow_lat_r + 0.5*elbow_med_r;

% wrist joint centers
wjc_l = 0.5*wrist_radial_l + 0.5*wrist_ulnar_l;
wjc_r = 0.5*wrist_radial_r + 0.5*wrist_ulnar_r;


%----------------------------------------------------------------
% Define the segments
%----------------------------------------------------------------

% The root segment: pelvis
pelvis.name = 'pelvis';
                                             
e_x = e_LR;
e_y = e_PA;
e_z = e_IS;

pelvis.localframe = cat(1, cat(2, e_x, e_y, e_z, midpelvis),...
			[0 0 0 1]);
pelvis.dof = {[1 2 3], [1 2 3]};

%states with typical range of motion
pelvis.states = {'pelvis x', 0.0821
		 'pelvis y', 0.0640
		 'pelvis z', 0.0770
		 'pelvis tilt', 0.1648
		 'pelvis obliqueity', 0.2129
		 'pelvis rotation', 0.3120};

% Tracking markers
pelvis.markers = {'Pelvis_1' pelvis_1
		  'Pelvis_2' pelvis_2
		  'Pelvis_3' pelvis_3};

%% Using LPT (lower part of trunk) data from de Leva
pelvis.length = 146*1e-3;
pelvis.mass = 0.112*mass;
pelvis.CoM = midpelvis;
pelvis.g0 = pelvis.localframe; %% Local coordinate system with origin at CoM of segment.
pelvis.moment_of_inertia = pelvis.mass ...
			   * diag( (pelvis.length*[0.615 0.551 0.587]).^2 );
pelvis.generalized_inertia = [pelvis.mass*eye(3) zeros(3,3)
			      zeros(3,3)  pelvis.moment_of_inertia];

% The trunk
trunk.name = 'trunk';
trunk_center = midpelvis; % Assume rotations around a point in
                              % the middle of the asis-psis plane.
trunk.localframe = cat(1, cat(2, e_x, e_y, e_z, trunk_center),...
		       [0 0 0 1]);
trunk.dof = {[2 1 3], []}; % The order of (euler) angles is y-x-z
trunk.states = {'trunk tilt', 0.2541
	        'trunk obliquety', 0.2468
	        'trunk rotation', 0.3226};
trunk.markers = {'Upper_Torso_1', ut_1
		 'Upper_Torso_2', ut_2
		 'Upper_Torso_3', ut_3};

%% Using MPT together with UPT (middle and upper part of trunk) data from de Leva
%% Ignoring the head, since it is close to still during the movement.
mpt.length = 216*1e-3;
mpt.mass = 0.163*mass;
mpt.CoM = midpelvis + 0.4*e_IS*mpt.length;
mpt.moment_of_inertia = mpt.mass  ...
			* diag( (mpt.length*[0.482 0.383 0.468]).^2 );
upt.length = 170*1e-3;
upt.mass = 0.16*mass;
upt.CoM = midtrunk - 0.3*e_IS*upt.length;
upt.moment_of_inertia = upt.mass  ...
			* diag( (upt.length*[0.716 0.454 0.659]).^2 );
trunk.mass = mpt.mass + upt.mass;
trunk.CoM = (mpt.CoM*mpt.mass + upt.CoM*upt.mass) / trunk.mass;
trunk.g0 = cat(1, cat(2, e_LR, e_PA, e_IS, trunk.CoM), [0 0 0 1]);
vm = mpt.CoM - trunk.CoM;
vu = upt.CoM - trunk.CoM;
trunk.moment_of_inertia = mpt.moment_of_inertia + mpt.mass * diag( [vm(2:3)'*vm(2:3)
								    vm([1 3])'*vm([1 3])
								    vm([1 2])'*vm([1 2])] ) ...
			  +upt.moment_of_inertia + upt.mass * diag( [vu(2:3)'*vu(2:3)
								    vu([1 3])'*vu([1 3])
								    vu([1 2])'*vu([1 2])] );
trunk.generalized_inertia = [trunk.mass*eye(3) zeros(3,3)
			     zeros(3,3) trunk.moment_of_inertia];

% The left arm 
luarm.name = 'left_upper_arm';
e_z = ejc_l - ghjc_l; % Local z-axis pointing axially from shoulder
                      % joint to elbow joint
e_z = e_z / norm(e_z);
e_x =  elbow_med_l - elbow_lat_l; % Local x-axis pointing
                                  % left-right 
e_x = e_x - (e_x'*e_z)*e_z;
e_x = e_x / norm(e_x);
e_y = cross(e_z, e_x);

luarm.localframe = cat(1, cat(2, e_x, e_y, e_z, ghjc_l),...
		      [0 0 0 1]);
luarm.dof = {[1 2 3], [1 2 3]};
luarm.states = {'left shoulder x', 0.1050
            'left shoulder y', 0.1148
            'left shoulder z', 0.0850
            'left shoulder flexion', 0.6668
            'left shoulder abduction', 0.7578
            'left shoulder rotation', 0.5297};
luarm.markers = {'L_Upper_Arm_1' , uarm_1_l
		'L_Upper_Arm_2' , uarm_2_l
		 'L_Upper_Arm_3' , uarm_3_l};

%% Data from de Leva
luarm.length = norm(ejc_l - ghjc_l);
luarm.mass = 0.027*mass;
luarm.CoM = ghjc_l + 0.577*luarm.length*e_z;
luarm.g0 = cat(1, cat(2, e_x, e_y, e_z, luarm.CoM), [0 0 0 1]);
luarm.moment_of_inertia = luarm.mass ...
			   * diag( (luarm.length*[0.285 0.269 0.158]).^2 );
luarm.generalized_inertia = [luarm.mass*eye(3) zeros(3,3)
			      zeros(3,3)  luarm.moment_of_inertia];

% Forearm
e_z = wjc_l - ejc_l;
e_z = e_z / norm(e_z);
e_x = e_x - (e_x'*e_z)*e_z; % Local x-axis similar to upper arm,
                            % but adjusted for the different axial direction
e_x = e_x / norm(e_x);
e_y = cross(e_z, e_x);

llarm.name = 'left_forearm';
llarm.localframe = cat(1, cat(2, e_x, e_y, e_z, ejc_l),...
		      [0 0 0 1]);
llarm.dof = {[1 3], []};
llarm.states = {'left elbow flexion', 0.5708
	      'left elbow rotation', 0.6794};
llarm.markers = {}; % No markers to track the forearm

%% Data from de Leva
llarm.length = norm(ejc_l - wjc_l);
llarm.mass = 0.0162*mass;
llarm.CoM = ejc_l + 0.457*llarm.length*e_z;
llarm.g0 = cat(1, cat(2, e_x, e_y, e_z, llarm.CoM), [0 0 0 1]);
llarm.moment_of_inertia = llarm.mass ...
			   * diag( (llarm.length*[0.276 0.265 0.121]).^2 );
llarm.generalized_inertia = [llarm.mass*eye(3) zeros(3,3)
			      zeros(3,3)  llarm.moment_of_inertia];

% Hand
mid_mp_l = 0.5*mp_5_l + 0.5*mp_2_l;

e_z = mid_mp_l - wjc_l;
e_z = e_z / norm(e_z);
e_x = mp_2_l - mp_5_l; % Local x-axis pointing left-right given by
                       % MP markers
e_x = e_x / norm(e_x);
% For the hand, adjust z-direction, not x-axis.
e_z = e_z - (e_z'*e_x)*e_x; 
e_z = e_z / norm(e_z);

e_y = cross(e_z, e_x);

lhand.name = 'left_hand';
lhand.localframe = cat(1, cat(2, e_x, e_y, e_z, wjc_l),...
		      [0 0 0 1]);
lhand.dof = {[1 2], []};
lhand.states = {'left wrist flexion', 0.5369
	      'left wrist abduction', 0.8487};
lhand.markers = {...
		'L_Hand_1', hand_1_l
		'L_Hand_2', hand_2_l
		'L_Hand_3', hand_3_l}; 
%% Data from de Leva
lhand.length = norm(mid_mp_l - wjc_l);
lhand.mass = 0.0061*mass;
lhand.CoM = wjc_l + 0.79*lhand.length*e_z;
lhand.g0 = cat(1, cat(2, e_x, e_y, e_z, lhand.CoM), [0 0 0 1]);
lhand.moment_of_inertia = lhand.mass ...
			   * diag( (lhand.length*[0.628 0.513 0.40]).^2 );
lhand.generalized_inertia = [lhand.mass*eye(3) zeros(3,3)
			      zeros(3,3)  lhand.moment_of_inertia];


% The right arm 
ruarm.name = 'right_upper_arm';
e_z = ejc_r - ghjc_r; % Local z-axis pointing axially from shoulder
                      % joint to elbow joint
e_z = e_z / norm(e_z);
e_x =  elbow_lat_r - elbow_med_r; % Local x-axis pointing
                                  % left-right 
e_x = e_x - (e_x'*e_z)*e_z;
e_x = e_x / norm(e_x);
e_y = cross(e_z, e_x);

ruarm.localframe = cat(1, cat(2, e_x, e_y, e_z, ghjc_r),...
		      [0 0 0 1]);
ruarm.dof = {[1 2 3], [1 2 3]};
ruarm.states = {'right shoulder x', 0.2678
            'right shoulder y', 0.1001
            'right shoulder z', 0.1048
            'right shoulder flexion', 0.5174
	      'right shoulder abduction', 0.3660
	       'right shoulder rotation', 1};
ruarm.markers = {'R_Upper_Arm_1' , uarm_1_r
		'R_Upper_Arm_2' , uarm_2_r
		 'R_Upper_Arm_3' , uarm_3_r};

%% Data from de Leva
ruarm.length = norm(ejc_r - ghjc_r);
ruarm.mass = 0.027*mass;
ruarm.CoM = ghjc_r + 0.577*ruarm.length*e_z;
ruarm.g0 = cat(1, cat(2, e_x, e_y, e_z, ruarm.CoM), [0 0 0 1]);
ruarm.moment_of_inertia = ruarm.mass ...
			   * diag( (ruarm.length*[0.285 0.269 0.158]).^2 );
ruarm.generalized_inertia = [ruarm.mass*eye(3) zeros(3,3)
			      zeros(3,3)  ruarm.moment_of_inertia];


% Forearm
e_z = wjc_r - ejc_r;
e_z = e_z / norm(e_z);
e_x = e_x - (e_x'*e_z)*e_z; % Local x-axis similar to upper arm,
                            % but adjusted for the different axial direction
e_x = e_x / norm(e_x);
e_y = cross(e_z, e_x);

rlarm.name = 'right_forearm';
rlarm.localframe = cat(1, cat(2, e_x, e_y, e_z, ejc_r),...
		      [0 0 0 1]);
rlarm.dof = {[1 3], []};
rlarm.states = {'right elbow flexion', 0.7490
	      'right elbow rotation', 0.2101};
rlarm.markers = {}; % No markers to track the forearm

%% Data from de Leva
rlarm.length = norm(ejc_r - wjc_r);
rlarm.mass = 0.0162*mass;
rlarm.CoM = ejc_r + 0.457*rlarm.length*e_z;
rlarm.g0 = cat(1, cat(2, e_x, e_y, e_z, rlarm.CoM), [0 0 0 1]);
rlarm.moment_of_inertia = rlarm.mass ...
			   * diag( (rlarm.length*[0.276 0.265 0.121]).^2 );
rlarm.generalized_inertia = [rlarm.mass*eye(3) zeros(3,3)
			      zeros(3,3)  rlarm.moment_of_inertia];

% Hand
mid_mp_r = 0.5*mp_5_r + 0.5*mp_2_r;

e_z = mid_mp_r - wjc_r;
e_z = e_z / norm(e_z);
e_x = mp_5_r - mp_2_r; % Local x-axis pointing left-right given by
                       % MP markers
e_x = e_x / norm(e_x);
% For the hand, adjust z-direction, not x-axis.
e_z = e_z - (e_z'*e_x)*e_x; 
e_z = e_z / norm(e_z);

e_y = cross(e_z, e_x);

rhand.name = 'right_hand';
rhand.localframe = cat(1, cat(2, e_x, e_y, e_z, wjc_r),...
		      [0 0 0 1]);
rhand.dof = {[1 2], []};
rhand.states = {'right wrist flexion', 0.9170
	      'right wrist abduction', 0.5507};
rhand.markers = {...
		'R_Hand_1', hand_1_r
		'R_Hand_2', hand_2_r
		'R_Hand_3', hand_3_r}; 

%% Data from de Leva
rhand.length = norm(mid_mp_r - wjc_r);
rhand.mass = 0.0061*mass;
rhand.CoM = wjc_r + 0.79*rhand.length*e_z;
rhand.g0 = cat(1, cat(2, e_x, e_y, e_z, rhand.CoM), [0 0 0 1]);
rhand.moment_of_inertia = rhand.mass ...
			   * diag( (rhand.length*[0.628 0.513 0.40]).^2 );
rhand.generalized_inertia = [rhand.mass*eye(3) zeros(3,3)
			      zeros(3,3)  rhand.moment_of_inertia];


% The club

% We need to define two club segments, but tracked using the same
% markers. Each "club" is at the end of either left and right
% kinematic chain.

e_z = grip_top_l - heel_bottom_grove_l;
e_z = e_z / norm(e_z);
e_y = toe_bottom_grove_l - heel_bottom_grove_l;
e_y = e_y - (e_y'*e_z)*e_z; 
e_y = e_y / norm(e_y);
e_x = cross(e_y, e_z);


% The club origin is taken as the point on the shaft closest to the
% wrist joint center
grip2hand = wjc_l - grip_top_l;
% Projection of the above vector onto shaft is the desired point
club_center = grip_top_l + (grip2hand'*e_z)*e_z;

club_l.name = 'club_l';
club_l.localframe = cat(1, cat(2, e_x, e_y, e_z, club_center),...
		      [0 0 0 1]);
club_l.dof = {[1 2 3], [1 2 3]};
club_l.states = {'club_l x', 0.1790
	       'club_l y', 0.1766
	       'club_l z', 0.0390
	       'club_l tilt', 0.5350
	       'club_l yaw', 0.6230
	       'club_l rotation', 0.2597};
club_l.markers = {...
		'Club_1', club_1_l
		'Club_2', club_2_l
		'Club_3', club_3_l}; 
club_l.CoM = 0.5*heel_bottom_grove_l + 0.5*toe_top_grove_l; % This point
                                                      % is used to
                                                      % compute the
                                                      % velocity of
                                                      % the club head
%% The inertial parameters. We divide the mass of the club in two for the "two" clubs.

%% Data from de Leva
rhand.length = norm(mid_mp_r - wjc_r);
rhand.mass = 0.0061*mass;
rhand.CoM = wjc_r + 0.79*rhand.length*e_z;
rhand.g0 = cat(1, cat(2, e_x, e_y, e_z, rhand.CoM), [0 0 0 1]);
rhand.moment_of_inertia = rhand.mass ...
			   * diag( (rhand.length*[0.628 0.513 0.40]).^2 );
rhand.generalized_inertia = [rhand.mass*eye(3) zeros(3,3)
			      zeros(3,3)  rhand.moment_of_inertia];



% The "right" club						      
e_z = grip_top_r - heel_bottom_grove_r;
e_z = e_z / norm(e_z);
e_y = toe_bottom_grove_r - heel_bottom_grove_r;
e_y = e_y - (e_y'*e_z)*e_z; 
e_y = e_y / norm(e_y);
e_x = cross(e_y, e_z);


% The club origin is taken as the point on the shaft closest to the
% wrist joint center
grip2hand = wjc_r - grip_top_r;
% Projection of the above vector onto shaft is the desired point
club_center = grip_top_r + (grip2hand'*e_z)*e_z;

club_r.name = 'club_r';
club_r.localframe = cat(1, cat(2, e_x, e_y, e_z, club_center),...
		      [0 0 0 1]);
club_r.dof = {[1 2 3], [1 2 3]};
club_r.states = {'club_r x', 0.1458
	       'club_r y', 0.1381
	       'club_r z', 0.0853
	       'club_r tilt', 0.3372
	       'club_r yaw', 0.3945
	       'club_r rotation', 0.1299};
club_r.markers = {...
		'Club_1', club_1_r
		'Club_2', club_2_r
		'Club_3', club_3_r}; 
club_r.CoM = 0.5*heel_bottom_grove_r + 0.5*toe_top_grove_r; % This point
                                                      % is used to
                                                      % compute the
                                                      % velocity of
                                                      % the club head
%club_r.CoM = club_1_r;
endpointstr = 'ClubCoM';


%----------------------------------------------------------------
% Define the complete model
%----------------------------------------------------------------

[tws, p0, gcnames, jc, segmnames, CoM] = build_model(pelvis);
[tws_ub, p0_ub, gcnames_ub, jc_ub, segmnames_ub, CoM_ub] = ...
    build_model(trunk);
[tws_la, p0_la, gcnames_la, jc_la, segmnames_la, CoM_la] = ...
    build_model(luarm, llarm, lhand, club_l);
[tws_ra, p0_ra, gcnames_ra, jc_ra, segmnames_ra, CoM_ra] = ...
    build_model(ruarm, rlarm, rhand, club_r);

tws_ub{2} = tws_la;
tws_ub{3} = tws_ra;

tws{2} = tws_ub;

p0_ub{2} = p0_la;
p0_ub{3} = p0_ra;

p0{2} = p0_ub;

jc_ub{2} = jc_la;
jc_ub{3} = jc_ra;

jc{2} = jc_ub;

gcnames_ub = cat(1, gcnames_ub, gcnames_la, gcnames_ra); 
gcnames = cat(1, gcnames, gcnames_ub);

segmnames_ub = cat(1, segmnames_ub, segmnames_la, segmnames_ra); 
segmnames = cat(1, segmnames, segmnames_ub);

CoM_ub{2} = CoM_la;
CoM_ub{3} = CoM_ra;

CoM{2} = CoM_ub;

tm.twists = tws;
tm.p0 = p0;
tm.jcs = jc;
tm.gcnames = gcnames;
tm.segm_names = segmnames;
tm.CoM = CoM;

function m = mmyextractmeanmarkers(rd, mname)
% Will look in struct rd for marker (or landmark) of name
% mname. Returns the average position in a 3x1 column vector, or
% NaNs if not found.

m = nan(3,1);

if isstruct(rd)
  
  if isfield(rd,mname)
    md = getfield(rd,mname);
    
    m = (mean(md{1}(:,1:3),1))';
  end

else
  mm = extractmarkers(rd,mname);
  m = mm(1,1:3)';
end
