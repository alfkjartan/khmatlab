function [gmleft, gmright, gmboth] = build_models_hand_as_endpoint(refdata, trialdata, bodymass)
%  [gmleft, gmright, gmboth] = build_models_hand_as_endpoint(refdata, trialdata, bodymass)
%
% Returns kinematic models for mobility calculations
%   gmleft     <-  Consists of the left arm with hand as endsegment.
%   gmright    <-  Consists of the right arm with hand as endsegment.
%   gmboth     <-  Both arms. 
%
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
if isempty(club_1_tr)
  club_1_tr = extractmarkers(trialdata, 'AnotherMarkerName');
end

[impact, impact_fit, pquad, dist2address, max_before_backsw] ...
= find_impact_from_point(club_1_tr);

trialdata{2} = trialdata{2}(max_before_backsw:max_before_backsw+1,:);
refdata = trialdata;

c7 = mmyextractmeanmarkers(refdata, 'C7');
ij = mmyextractmeanmarkers(refdata, 'Insicura Jugularis');
shoulder_l = mmyextractmeanmarkers(refdata, 'L Acromion');
shoulder_r = mmyextractmeanmarkers(refdata, 'R Acromion');
ghjc_l = mmyextractmeanmarkers(refdata, 'Wrt_LShoulder');
ghjc_r = mmyextractmeanmarkers(refdata, 'Wrt_RShoulder');
elbow_lat_l = mmyextractmeanmarkers(refdata, 'L Elbow lateral');
elbow_med_l = mmyextractmeanmarkers(refdata, 'L Elbow medial');
elbow_lat_r = mmyextractmeanmarkers(refdata, 'R Elbow lateral');
elbow_med_r = mmyextractmeanmarkers(refdata, 'R Elbow medial');
wrist_radial_l = mmyextractmeanmarkers(refdata, 'L Radial wrist');
wrist_ulnar_l = mmyextractmeanmarkers(refdata, 'L Ulnar wrist');
wrist_radial_r = mmyextractmeanmarkers(refdata, 'R Radial wrist');
wrist_ulnar_r = mmyextractmeanmarkers(refdata, 'R Ulnar wrist');
asis_l = mmyextractmeanmarkers(refdata, 'L ASIS');
asis_r = mmyextractmeanmarkers(refdata, 'R ASIS');
psis_l = mmyextractmeanmarkers(refdata, 'L PSIS');
psis_r = mmyextractmeanmarkers(refdata, 'R PSIS');
t8 = mmyextractmeanmarkers(refdata, 'T8');
pelvis_1 = mmyextractmeanmarkers(refdata, 'PELVIS_1');
pelvis_2 = mmyextractmeanmarkers(refdata, 'PELVIS_2');
pelvis_3 = mmyextractmeanmarkers(refdata, 'PELVIS_3');
ut_1 = mmyextractmeanmarkers(refdata, 'UPPER_TORSO_1');
ut_2 = mmyextractmeanmarkers(refdata, 'UPPER_TORSO_2');
ut_3 = mmyextractmeanmarkers(refdata, 'UPPER_TORSO_3');
uarm_1_l = mmyextractmeanmarkers(refdata, 'L_UPPER_ARM_1');
uarm_2_l = mmyextractmeanmarkers(refdata, 'L_UPPER_ARM_2');
uarm_3_l = mmyextractmeanmarkers(refdata, 'L_UPPER_ARM_3');
uarm_1_r = mmyextractmeanmarkers(refdata, 'R_UPPER_ARM_1');
uarm_2_r = mmyextractmeanmarkers(refdata, 'R_UPPER_ARM_2');
uarm_3_r = mmyextractmeanmarkers(refdata, 'R_UPPER_ARM_3');
hand_1_l = mmyextractmeanmarkers(refdata, 'L_HAND_1');
hand_2_l = mmyextractmeanmarkers(refdata, 'L_HAND_2');
hand_3_l = mmyextractmeanmarkers(refdata, 'L_HAND_3');
hand_1_r = mmyextractmeanmarkers(refdata, 'R_HAND_1');
hand_2_r = mmyextractmeanmarkers(refdata, 'R_HAND_2');
hand_3_r = mmyextractmeanmarkers(refdata, 'R_HAND_3');
hand_1_l_trial = mmyextractmeanmarkers(trialdata, 'L_HAND_1');
hand_2_l_trial = mmyextractmeanmarkers(trialdata, 'L_HAND_2');
hand_3_l_trial = mmyextractmeanmarkers(trialdata, 'L_HAND_3');
hand_1_r_trial = mmyextractmeanmarkers(trialdata, 'R_HAND_1');
hand_2_r_trial = mmyextractmeanmarkers(trialdata, 'R_HAND_2');
hand_3_r_trial = mmyextractmeanmarkers(trialdata, 'R_HAND_3');
mp_2_l = mmyextractmeanmarkers(refdata, 'L 2nd MP joint');
mp_5_l = mmyextractmeanmarkers(refdata, 'L 5th MP joint');
mp_2_r = mmyextractmeanmarkers(refdata, 'R 2nd MP joint');
mp_5_r = mmyextractmeanmarkers(refdata, 'R 5th MP joint');
grip_top = mmyextractmeanmarkers(refdata, 'Top of handle');
heel_bottom_grove = mmyextractmeanmarkers(refdata, 'Bottom grove@heel');
toe_bottom_grove = mmyextractmeanmarkers(refdata, 'Bottom grove@toe');
toe_top_grove = mmyextractmeanmarkers(refdata, 'Top grove@toe');
club_1 = mmyextractmeanmarkers(refdata, 'Club_1');
club_2 = mmyextractmeanmarkers(refdata, 'Club_2');
club_3 = mmyextractmeanmarkers(refdata, 'Club_3');
club_1_trial = mmyextractmeanmarkers(trialdata, 'CLUB_1');
club_2_trial = mmyextractmeanmarkers(trialdata, 'CLUB_2');
club_3_trial = mmyextractmeanmarkers(trialdata, 'CLUB_3');
%midhands_1 = mmyextractmeanmarkers(refdata, 'MidHands_1');
%midhands_2 = mmyextractmeanmarkers(refdata, 'MidHands_2');

		       
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
lhand.dof = {[1 2 ], []};
lhand.states = {'left wrist flexion', 0.9170
		'left wrist abduction', 0.5507};
lhand.markers = {'L_Hand_1', hand_1_l
		 'L_Hand_2', hand_2_l
		 'L_Hand_3', hand_3_l
		 'CLUB_1', club_1
		 'CLUB_2', club_2
		 'CLUB_3', club_3};
%% Data from de Leva
lhand.length = norm(mid_mp_l - wjc_l);
lhand.mass = 0.0061*bodymass;
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
		 'CLUB_1', club_1
		 'CLUB_2', club_2
		 'CLUB_3', club_3};
%% Data from de Leva
rhand.length = norm(mid_mp_r - wjc_r);
rhand.mass = 0.0061*bodymass;
rhand.CoM = wjc_r + 0.79*rhand.length*e_z;
rhand.g0 = cat(1, cat(2, e_x, e_y, e_z, rhand.CoM), [0 0 0 1]);
rhand.moment_of_inertia = rhand.mass ...
			   * diag( (rhand.length*[0.628 0.513 0.40]).^2 );
rhand.generalized_inertia = [rhand.mass*eye(3) zeros(3,3)
			      zeros(3,3)  rhand.moment_of_inertia];

% Object frame is club frame
e_z = grip_top - heel_bottom_grove;
e_z = e_z / norm(e_z);
e_y = toe_bottom_grove - heel_bottom_grove;
e_y = e_y - (e_y'*e_z)*e_z; 
e_y = e_y / norm(e_y);
e_x = cross(e_y, e_z);

endpoint =0.5*rhand.CoM + 0.5*lhand.CoM; 
% This point is used to compute the velocity of the club head
club_g0 = cat(1, cat(2, e_x, e_y, e_z, endpoint),...
			[0 0 0 1]);
lhand.object_frame = club_g0;
rhand.object_frame = club_g0;

%%keyboard
endpointstr = 'MidHands';


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
