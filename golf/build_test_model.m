function tm = build_test_model(refdata,trialdata)
%  gm = build_test_model(refdata, trialdata)
%
% Returns a kinematic model for testing the anglecontribution
% function. 
% The model 2 segments: pelvis and trunk
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
% Output
%    tm        <-  struct. Contains the fields 
%        tws       <-  nested cell array of twists
%        gcnames   <-  array with names for the generalized
%                      coordinates.
%        p0        <-  nested cell with reference marker positions.
%        jcs       <-  nested cell with joint centra 
%        CoM       <-  nested cell with center of mass 
%

% Kjartan Halvorsen
% 2009-06-26
% Based on build_tennis_model.m from  2003-09-16
  
c7 = myextractmeanmarkers(refdata, 'C7');
ij = myextractmeanmarkers(refdata, 'IJ');
shoulder_l = myextractmeanmarkers(refdata, 'L_Acromion');
shoulder_r = myextractmeanmarkers(refdata, 'R_Acromion');
ghjc_l = myextractmeanmarkers(refdata, 'Wrt_LShoulder');
ghjc_r = myextractmeanmarkers(refdata, 'Wrt_RShoulder');
elbow_lat_l = myextractmeanmarkers(refdata, 'L_Elbow_lateral');
elbow_med_l = myextractmeanmarkers(refdata, 'L_Elbow_medial');
elbow_lat_r = myextractmeanmarkers(refdata, 'R_Elbow_lateral');
elbow_med_r = myextractmeanmarkers(refdata, 'R_Elbow_medial');
wrist_radial_l = myextractmeanmarkers(refdata, 'L_Radial_wrist');
wrist_ulnar_l = myextractmeanmarkers(refdata, 'L_Ulnar_wrist');
wrist_radial_r = myextractmeanmarkers(refdata, 'R_Radial_wrist');
wrist_ulnar_r = myextractmeanmarkers(refdata, 'R_Ulnar_wrist');
asis_l = myextractmeanmarkers(refdata, 'L_ASIS');
asis_r = myextractmeanmarkers(refdata, 'R_ASIS');
psis_l = myextractmeanmarkers(refdata, 'L_PSIS');
psis_r = myextractmeanmarkers(refdata, 'R_PSIS');
t8 = myextractmeanmarkers(refdata, 'T8');
pelvis_1 = myextractmeanmarkers(refdata, 'Pelvis_1');
pelvis_2 = myextractmeanmarkers(refdata, 'Pelvis_2');
pelvis_3 = myextractmeanmarkers(refdata, 'Pelvis_3');
ut_1 = myextractmeanmarkers(refdata, 'Upper_Torso_1');
ut_2 = myextractmeanmarkers(refdata, 'Upper_Torso_2');
ut_3 = myextractmeanmarkers(refdata, 'Upper_Torso_3');
uarm_1_l = myextractmeanmarkers(refdata, 'L_Upper_Arm_1');
uarm_2_l = myextractmeanmarkers(refdata, 'L_Upper_Arm_2');
uarm_3_l = myextractmeanmarkers(refdata, 'L_Upper_Arm_3');
uarm_1_r = myextractmeanmarkers(refdata, 'R_Upper_Arm_1');
uarm_2_r = myextractmeanmarkers(refdata, 'R_Upper_Arm_2');
uarm_3_r = myextractmeanmarkers(refdata, 'R_Upper_Arm_3');
hand_1_l = myextractmeanmarkers(refdata, 'L_Hand_1');
hand_2_l = myextractmeanmarkers(refdata, 'L_Hand_2');
hand_3_l = myextractmeanmarkers(refdata, 'L_Hand_3');
hand_1_r = myextractmeanmarkers(refdata, 'R_Hand_1');
hand_2_r = myextractmeanmarkers(refdata, 'R_Hand_2');
hand_3_r = myextractmeanmarkers(refdata, 'R_Hand_3');
hand_1_l_trial = myextractmeanmarkers(trialdata, 'L_Hand_1');
hand_2_l_trial = myextractmeanmarkers(trialdata, 'L_Hand_2');
hand_3_l_trial = myextractmeanmarkers(trialdata, 'L_Hand_3');
hand_1_r_trial = myextractmeanmarkers(trialdata, 'R_Hand_1');
hand_2_r_trial = myextractmeanmarkers(trialdata, 'R_Hand_2');
hand_3_r_trial = myextractmeanmarkers(trialdata, 'R_Hand_3');
mp_2_l = myextractmeanmarkers(refdata, 'L_2nd_MP_joint');
mp_5_l = myextractmeanmarkers(refdata, 'L_5th_MP_joint');
mp_2_r = myextractmeanmarkers(refdata, 'R_2nd_MP_joint');
mp_5_r = myextractmeanmarkers(refdata, 'R_5th_MP_joint');
grip_top = myextractmeanmarkers(refdata, 'TopOfHandle');
heel_bottom_grove = myextractmeanmarkers(refdata, 'BottomGroveHeel');
toe_bottom_grove = myextractmeanmarkers(refdata, 'BottomGroveToe');
toe_top_grove = myextractmeanmarkers(refdata, 'TopGroveToe');
club_1 = myextractmeanmarkers(refdata, 'Club_1');
club_2 = myextractmeanmarkers(refdata, 'Club_2');
club_3 = myextractmeanmarkers(refdata, 'Club_3');
club_1_trial = myextractmeanmarkers(trialdata, 'Club_1');
club_2_trial = myextractmeanmarkers(trialdata, 'Club_2');
club_3_trial = myextractmeanmarkers(trialdata, 'Club_3');
%midhands_1 = myextractmeanmarkers(refdata, 'MidHands_1');
%midhands_2 = myextractmeanmarkers(refdata, 'MidHands_2');

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
pelvis.states = {'pelvis x', 0.2
		 'pelvis y', 0.2
		 'pelvis z', 0.1
		 'pelvis tilt', pi/3
		 'pelvis obliqueity', pi/3
		 'pelvis rotation', pi};

% Tracking markers
pelvis.markers = {'Pelvis_1' pelvis_1
		  'Pelvis_2' pelvis_2
		  'Pelvis_3' pelvis_3};

% pelvis.mass = 0.112;


% The trunk
trunk.name = 'trunk';
trunk_center = midpelvis; % Assume rotations around a point in
                              % the middle of the asis-psis plane.
trunk.localframe = cat(1, cat(2, e_x, e_y, e_z, trunk_center),...
		       [0 0 0 1]);
trunk.dof = {[2 1 3], []}; % The order of (euler) angles is y-x-z
trunk.states = {'trunk tilt', pi/2
	        'trunk obliquety', pi/6
	        'trunk rotation', pi/2};
trunk.markers = {'Upper_Torso_1', ut_1
		 'Upper_Torso_2', ut_2
		 'Upper_Torso_3', ut_3};

trunk.CoM = ut_1;
                                                      % is used to
                                                      % compute the
                                                      % velocity of
                                                      % the club head

if 0
  %trunk.mass = 0.323;

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
  luarm.dof = {[1 2 3], []};
  luarm.states = {'left shoulder flexion', pi/2
		  'left shoulder abduction', pi/2
		  'left shoulder rotation', pi/2};
  luarm.markers = {'L_Upper_Arm_1' , uarm_1_l
		   'L_Upper_Arm_2' , uarm_2_l
		   'L_Upper_Arm_3' , uarm_3_l};

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
  llarm.states = {'left elbow flexion', pi/2
		  'left elbow rotation', pi/2};
  llarm.markers = {}; % No markers to track the forearm

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
  lhand.states = {'left wrist flexion', pi/2
		  'left wrist abduction', pi/2};
  lhand.markers = {...
      'L_Hand_1', hand_1_l
      'L_Hand_2', hand_2_l
      'L_Hand_3', hand_3_l}; 


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
  ruarm.dof = {[1 2 3], []};
  ruarm.states = {'right shoulder flexion', pi/2
		  'right shoulder abduction', pi/2
		  'right shoulder rotation', pi/2};
  ruarm.markers = {'R_Upper_Arm_1' , uarm_1_r
		   'R_Upper_Arm_2' , uarm_2_r
		   'R_Upper_Arm_3' , uarm_3_r};

  % Forearm
  e_z = wjc_l - ejc_l;
  e_z = e_z / norm(e_z);
  e_x = e_x - (e_x'*e_z)*e_z; % Local x-axis similar to upper arm,
                            % but adjusted for the different axial direction
  e_x = e_x / norm(e_x);
  e_y = cross(e_z, e_x);

  rlarm.name = 'right_forearm';
  rlarm.localframe = cat(1, cat(2, e_x, e_y, e_z, ejc_r),...
			 [0 0 0 1]);
  rlarm.dof = {[1 3], []};
  rlarm.states = {'right elbow flexion', pi/2
		  'right elbow rotation', pi/2};
  rlarm.markers = {}; % No markers to track the forearm

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
  rhand.states = {'right wrist flexion', pi/2
		  'right wrist abduction', pi/2};
  rhand.markers = {...
      'R_Hand_1', hand_1_r
      'R_Hand_2', hand_2_r
      'R_Hand_3', hand_3_r}; 
  
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
  club_l.states = {'club_l x', 0.1
		   'club_l y', 0.1
		   'club_l z', 0.1
		   'club_l tilt', pi/2
		   'club_l yaw', pi/2 
		   'club_l rotation', pi/2};
  club_l.markers = {...
      'Club_1', club_1_l
      'Club_2', club_2_l
      'Club_3', club_3_l}; 
  club_l.CoM = 0.5*heel_bottom_grove_l + 0.5*toe_top_grove_l; % This point
                                                      % is used to
                                                      % compute the
                                                      % velocity of
                                                      % the club head

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
  club_r.states = {'club_r x', 0.1
		   'club_r y', 0.1
		   'club_r z', 0.1
		   'club_r tilt', pi/2
		   'club_r yaw', pi/2 
		   'club_r rotation', pi/2};
  club_r.markers = {...
      'Club_1', club_1_r
      'Club_2', club_2_r
      'Club_3', club_3_r}; 
  club_r.CoM = 0.5*heel_bottom_grove_r + 0.5*toe_top_grove_r; % This point
							      % is used to
							      % compute the
							      % velocity of
							      % the club head
end % if 0

%----------------------------------------------------------------
% Define the complete model
%----------------------------------------------------------------

[tws, p0, gcnames, jc, segmnames, CoM] = build_model(pelvis);
[tws_ub, p0_ub, gcnames_ub, jc_ub, segmnames_ub, CoM_ub] = ...
    build_model(trunk);
%[tws_la, p0_la, gcnames_la, jc_la, segmnames_la, CoM_la] = ...
%    build_model(luarm, llarm, lhand, club_l);
%[tws_ra, p0_ra, gcnames_ra, jc_ra, segmnames_ra, CoM_ra] = ...
%    build_model(ruarm, rlarm, rhand, club_r);

%tws_ub{2} = tws_la;
%tws_ub{3} = tws_ra;

tws{2} = tws_ub;

%p0_ub{2} = p0_la;
%p0_ub{3} = p0_ra;

p0{2} = p0_ub;

%jc_ub{2} = jc_la;
%jc_ub{3} = jc_ra;

jc{2} = jc_ub;

%gcnames_ub = cat(1, gcnames_ub, gcnames_la, gcnames_ra); 
gcnames = cat(1, gcnames, gcnames_ub);

%segmnames_ub = cat(1, segmnames_ub, segmnames_la, segmnames_ra); 
segmnames = cat(1, segmnames, segmnames_ub);

%CoM_ub{2} = CoM_la;
%CoM_ub{3} = CoM_ra;

CoM{2} = CoM_ub;

tm.twists = tws;
tm.p0 = p0;
tm.jcs = jc;
tm.gcnames = gcnames;
tm.segm_names = segmnames;
tm.CoM = CoM;

function m = myextractmeanmarkers(rd, mname)
% Will look in struct rd for marker (or landmark) of name
% mname. Returns the average position in a 3x1 column vector, or
% NaNs if not found.

m = nan(3,1);

if isstruct(rd)
  
  if isfield(rd,mname)
    md = getfield(rd,mname);
    m = (mean(md{1}(:,1:3)))';
  end

else
  mm = extractmarkers(rd,mname);
  m = mm(1,1:3)';
end
