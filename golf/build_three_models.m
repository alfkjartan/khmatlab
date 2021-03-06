function [gmleft, gmright, gmboth, gmclub] = build_three_models(refdata, trialdata, bodymass)
%  [gmleft, gmright, gmboth, gmclub] = build_three_models(refdata, trialdata, bodymass)
%
% Returns three kinematic models for the golf study:
%   gmleft     <-  Consists of the left upper arm, forearm and hand
%   gmright    <-  Consists of the right upper arm, forearm and hand
%   gmboth     <-  The two above, without the club. The inertia of the club needs to be
%                  taken into account separately.
%   gmclub     <-  Model of the club.
% Note that the grip is assumed to be firm, so that the club and hand(s) form one segment.
% In other words: no degrees of freedom between the hands and club.
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
% Output
%    gm{left,right,both}  <-  struct. Contains the fields 
%        tws       <-  nested cell array of twists
%        gcnames   <-  array with names for the generalized
%                      coordinates.
%        p0        <-  nested cell array with reference marker positions.
%        jcs       <-  nested cell array with joint centra 
%        CoM       <-  nested cell array with center of mass 
%        inertia   <-  nested cell array with inertial matrices (6 x 6)
%        g0        <-  nested cell array with local coordinate frames (6 x 6)
%        object_frame  <-  nested cell array with local object frames (6 x 6). 
%                          Typically only given for the end segment.

% Kjartan Halvorsen
% 2013-08-23
% Based on build_golf_model_w_inertia
  
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
hand_1_l_trial = mmyextractmeanmarkers(trialdata, 'L_HAND_1');
hand_2_l_trial = mmyextractmeanmarkers(trialdata, 'L_HAND_2');
hand_3_l_trial = mmyextractmeanmarkers(trialdata, 'L_HAND_3');
hand_1_r_trial = mmyextractmeanmarkers(trialdata, 'R_HAND_1');
hand_2_r_trial = mmyextractmeanmarkers(trialdata, 'R_HAND_2');
hand_3_r_trial = mmyextractmeanmarkers(trialdata, 'R_HAND_3');
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
club_1_trial = mmyextractmeanmarkers(trialdata, 'CLUB_1');
club_2_trial = mmyextractmeanmarkers(trialdata, 'CLUB_2');
club_3_trial = mmyextractmeanmarkers(trialdata, 'CLUB_3');
%midhands_1 = mmyextractmeanmarkers(refdata, 'MidHands_1');
%midhands_2 = mmyextractmeanmarkers(refdata, 'MidHands_2');



% Transform the club markers so that they appear in the reference trial (which is
% as always neutral standing) as if gripped as in addressing the ball. This means that
% there will be separate reference club points for each arm. 
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
%            'left shoulder flexion', 0.6668
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
luarm.mass = 0.027*bodymass;
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
llarm.mass = 0.0162*bodymass;
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
lhand.markers = {'L_Hand_1', hand_1_l
		 'L_Hand_2', hand_2_l
		 'L_Hand_3', hand_3_l
		 'CLUB_1', club_1_l
		 'CLUB_2', club_2_l
		 'CLUB_3', club_3_l};
%% Data from de Leva
lhand.length = norm(mid_mp_l - wjc_l);
lhand.mass = 0.0061*bodymass;
lhand.CoM = wjc_l + 0.79*lhand.length*e_z;
lhand.g0 = cat(1, cat(2, e_x, e_y, e_z, lhand.CoM), [0 0 0 1]);
lhand.moment_of_inertia = lhand.mass ...
			   * diag( (lhand.length*[0.628 0.513 0.40]).^2 );
lhand.generalized_inertia = [lhand.mass*eye(3) zeros(3,3)
			      zeros(3,3)  lhand.moment_of_inertia];

% Object frame is left club frame
e_z = grip_top_l - heel_bottom_grove_l;
e_z = e_z / norm(e_z);
e_y = toe_bottom_grove_l - heel_bottom_grove_l;
e_y = e_y - (e_y'*e_z)*e_z; 
e_y = e_y / norm(e_y);
e_x = cross(e_y, e_z);

club_head_center = 0.5*heel_bottom_grove_l + 0.5*toe_top_grove_l; 
% This point is used to compute the velocity of the club head
lhand.object_frame = cat(1, cat(2, e_x, e_y, e_z, club_head_center),...
			[0 0 0 1]);



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
ruarm.mass = 0.027*bodymass;
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
rlarm.mass = 0.0162*bodymass;
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
rhand.markers = {'R_Hand_1', hand_1_r
		 'R_Hand_2', hand_2_r
		 'R_Hand_3', hand_3_r
		 'CLUB_1', club_1_r
		 'CLUB_2', club_2_r
		 'CLUB_3', club_3_r};
%% Data from de Leva
rhand.length = norm(mid_mp_r - wjc_r);
rhand.mass = 0.0061*bodymass;
rhand.CoM = wjc_r + 0.79*rhand.length*e_z;
rhand.g0 = cat(1, cat(2, e_x, e_y, e_z, rhand.CoM), [0 0 0 1]);
rhand.moment_of_inertia = rhand.mass ...
			   * diag( (rhand.length*[0.628 0.513 0.40]).^2 );
rhand.generalized_inertia = [rhand.mass*eye(3) zeros(3,3)
			      zeros(3,3)  rhand.moment_of_inertia];

% Object frame is right club frame
e_z = grip_top_r - heel_bottom_grove_r;
e_z = e_z / norm(e_z);
e_y = toe_bottom_grove_r - heel_bottom_grove_r;
e_y = e_y - (e_y'*e_z)*e_z; 
e_y = e_y / norm(e_y);
e_x = cross(e_y, e_z);

club_head_center = 0.5*heel_bottom_grove_r + 0.5*toe_top_grove_r; 
% This point is used to compute the velocity of the club head
rhand.object_frame = cat(1, cat(2, e_x, e_y, e_z, club_head_center),...
			[0 0 0 1]);



%% The club

% The club is assumed to be made up of three parts: grip, shaft and head. 
% The moment of inertia is determined for each of these parts, and combined
% to form a hand+club segment. The steps are:
%   1) Determine inertia matrix in local coordinate system of part.
%   2) Transfrom to local coordinate system of hand
%   3) Add to inertia of hand using the paralell axis theorem.
e_shaft = heel_bottom_grove - grip_top; % Points down the shaft
e_shaft = e_shaft / norm(e_shaft);
e_club_pa = toe_bottom_grove - heel_bottom_grove;
e_club_pa = e_club_pa - (e_club_pa'*e_shaft)*e_shaft;
e_club_pa = e_club_pa / norm(e_club_pa);
e_club_x = cross(e_club_pa, e_shaft);
R_club = [e_club_x e_club_pa e_shaft];

club_head_center = 0.5*heel_bottom_grove + 0.5*toe_top_grove; 
% This point is used to compute the velocity of the club head
club.object_frame = cat(1, cat(2, e_club_x, e_club_pa, e_shaft, club_head_center),...
			[0 0 0 1]);
grip.mass = 50*1e-3; % 50 g
grip_length = 0.27; % 27cm typical length of the grip
grip_r2 = 0.01; % Radius of outer contour of the grip
grip_r1 = 0.006; % 6mm radius of shaft 
grip.CoM = grip_top + 0.4*grip_length*e_shaft; % CoM 12cm down the shaft. Grip typically 27cm long.
grip.local_inertia = grip.mass * diag( [1/12*(3*(grip_r1^2 + grip_r2^2) + grip_length*2)
					1/12*(3*(grip_r1^2 + grip_r2^2) + grip_length*2)
					1/2 * (grip_r1^2 + grip_r2^2)]);
grip.inertia = R_club * grip.local_inertia * R_club';

shaft.mass = 120*1e-3; % 120 g
shaft_length = norm(grip_top_l - heel_bottom_grove_l); % length of the shaft
shaft_r = 0.004; % average radius of  shaft
shaft.CoM = grip_top + 0.35*shaft_length*e_shaft; 
shaft.local_inertia = shaft.mass * diag( [1/12*(6*shaft_r^2 + shaft_length*2)
					  1/12*(6*shaft_r^2 + shaft_length*2)
					  shaft_r^2]);
shaft.inertia = R_club * shaft.local_inertia * R_club';

%% Clubhead as point mass
clubhead.mass = 300*1e-3; % 300 g
clubhead.CoM = (0.5* toe_bottom_grove + 0.5*heel_bottom_grove);
clubhead.inertia = zeros(3,3);

%% Combine
[club.inertia, club.CoM, club.mass] = combine_inertia(grip.inertia, grip.CoM, grip.mass, shaft.inertia, shaft.CoM, shaft.mass);
[club.inertia, club.CoM, club.mass] = combine_inertia(club.inertia, club.CoM, club.mass, clubhead.inertia, clubhead.CoM, clubhead.mass);

club.moment_of_inertia = club.inertia;
club.generalized_inertia = [club.mass*eye(3) zeros(3,3)
			    zeros(3,3) club.moment_of_inertia];

club.g0 = cat(1, cat(2, e_club_x, e_club_pa, e_shaft, club.CoM),...
			[0 0 0 1]);
club.name = 'club';
club.localframe = club.object_frame;
club.dof = {[1 2 3], [1 2 3]};
club.states = { 'club ml trl', 1
	       'club ap trl', 0.5
	       'club ax trl', 0.5
	       'club ml rot', pi
	       'club ap rot', pi
	       'club ax rot', pi};

club.markers = { 'CLUB_1', club_1
		 'CLUB_2', club_2
		 'CLUB_3', club_3};

endpointstr = 'ClubCoM';


%%%%%%%%%%%%%%%%%%%%%
% Need also a dummy trunk segment to act as root
%%%%%%%%%%%%%%%%%%%%%

trunk.CoM = zeros(3,1);
trunk.mass = 0;
trunk.localframe = eye(4);
trunk.dof = {[], []};
trunk.states = {};
trunk.moment_of_inertia = zeros(3,3);
trunk.generalized_inertia = zeros(6,6);
trunk.g0 = eye(4);


%----------------------------------------------------------------
% Define the complete models
%----------------------------------------------------------------

[gmboth.twists, gmboth.p0, gmboth.gcnames, gmboth.jcs, gmboth.segm_names, gmboth.CoM, radius, gmboth.mass, gmboth.g0, gmboth.inertia, gmboth.object_frame, gmboth.objectcenter] = build_model(trunk);
[gmleft.twists, gmleft.p0, gmleft.gcnames, gmleft.jcs, gmleft.segm_names, gmleft.CoM, ...
 radius_la, gmleft.mass, gmleft.g0, gmleft.inertia, gmleft.object_frame, gmleft.objectcenter] =  build_model(luarm, llarm, lhand);
[gmright.twists, gmright.p0, gmright.gcnames, gmright.jcs, gmright.segm_names, gmright.CoM, ...
 radius_la, gmright.mass, gmright.g0, gmright.inertia, gmright.object_frame, gmright.objectcenter] =  build_model(ruarm, rlarm, rhand);

gmboth.twists{2} = gmleft.twists;
gmboth.twists{3} = gmright.twists;
gmboth.p0{2} = gmleft.p0;
gmboth.p0{3} = gmright.p0;
gmboth.jcs{2} = gmleft.jcs;
gmboth.jcs{3} = gmright.jcs;
gmboth.CoM{2} = gmleft.CoM;
gmboth.CoM{3} = gmright.CoM;
gmboth.g0{2} = gmleft.g0;
gmboth.g0{3} = gmright.g0;
gmboth.inertia{2} = gmleft.inertia;
gmboth.inertia{3} = gmright.inertia;
gmboth.object_frame{2} = gmleft.object_frame;
gmboth.object_frame{3} = gmright.object_frame;
gmboth.objectcenter{2} = gmleft.objectcenter;
gmboth.objectcenter{3} = gmright.objectcenter;
gmboth.gcnames = cat(1, gmboth.gcnames, gmleft.gcnames, gmright.gcnames);
gmboth.segm_names = cat(1, gmboth.segm_names, gmleft.segm_names, gmright.segm_names);


[gmclub.twists, gmclub.p0, gmclub.gcnames, gmclub.jcs, gmclub.segm_names, gmclub.CoM, radius, gmclub.mass, gmclub.g0, gmclub.inertia, gmclub.object_frame, gmclub.objectcenter] = build_model(club);

%keyboard

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
