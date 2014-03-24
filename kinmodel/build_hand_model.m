function [tm, endpointstr] = build_hand_model(refdata)
%  gm = build_golf_model(refdata)
%
% Returns a kinematic model for the hand tracking data
% The model contains 9 segments:
%    pelvis, trunk, upperarms, underarms, hands, club
% The pelvis is the root segment.
%
% Input
%    refata    ->  Reference marker data. 
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
% 2011-05-19
% Based on build_golf_model from  2009-06-26
  
%% Revisions
%% 2011-09-12  Changed CoM to middle of distal tip. Used for spherical
%%             model of tip to detect contact. Added field "radius" to
%%             hold radius of distal tip.

styl_rad = myextractmeanmarkers(refdata, 'styl_radii');
styl_uln = myextractmeanmarkers(refdata, 'styl_ulnae');
radius = myextractmeanmarkers(refdata, 'radius');
mc1 = myextractmeanmarkers(refdata, 'mc1');
mcp1 = myextractmeanmarkers(refdata, 'mcp1');
pip1 = myextractmeanmarkers(refdata, 'pip1');
tip1 = myextractmeanmarkers(refdata, 'tip1');
mcp2 = myextractmeanmarkers(refdata, 'mcp2');
pip2 = myextractmeanmarkers(refdata, 'pip2');
dip2 = myextractmeanmarkers(refdata, 'dip2');
tip2 = myextractmeanmarkers(refdata, 'tip2');
mcp3 = myextractmeanmarkers(refdata, 'mcp3');
pip3 = myextractmeanmarkers(refdata, 'pip3');
dip3 = myextractmeanmarkers(refdata, 'dip3');
tip3 = myextractmeanmarkers(refdata, 'tip3');
mcp4 = myextractmeanmarkers(refdata, 'mcp4');
pip4 = myextractmeanmarkers(refdata, 'pip4');
dip4 = myextractmeanmarkers(refdata, 'dip4');
tip4 = myextractmeanmarkers(refdata, 'tip4');
mcp5 = myextractmeanmarkers(refdata, 'mcp5');
pip5 = myextractmeanmarkers(refdata, 'pip5');
dip5 = myextractmeanmarkers(refdata, 'dip5');
tip5 = myextractmeanmarkers(refdata, 'tip5');
plate1 = myextractmeanmarkers(refdata, 'plate1');
plate2 = myextractmeanmarkers(refdata, 'plate2');
plate3 = myextractmeanmarkers(refdata, 'plate3');
plate4 = myextractmeanmarkers(refdata, 'plate4');
plate5 = myextractmeanmarkers(refdata, 'plate5');

%keyboard

% Landmarks and anatomical directions
% Medio-lateralt
e_ml = styl_uln - styl_rad;
e_ml = e_ml / norm(e_ml);




%----------------------------------------------------------------
% Define the segments
%----------------------------------------------------------------

%-------------------
% Forearm
%-------------------

% Local coordinate system of the forearm
e_x = e_ml;
e_z = radius - styl_rad;
e_z = e_z - (e_z'*e_x)*e_x;
e_z = e_z / norm(e_z);
e_y = cross(e_z, e_x);

e_pd = e_y; % plantar-dorsal direction

% wrist joint center
wjc = 0.5*styl_uln + 0.5*styl_rad - 0.01*e_y; 
% 1cm below line between stylodeii. 

% The root segment: forearm
forearm.name = 'forearm';
                                             
forearm.localframe = cat(1, cat(2, e_x, e_y, e_z, wjc),...
			[0 0 0 1]);
forearm.dof = {[1 2 3], [1 2 3]};

%states with typical range of motion
forearm.states = {'forearm x', 0.2
		 'forearm y', 0.2
		 'forearm z', 0.1
		 'forearm flex', pi/3
		 'forearm add', pi/3
		 'forearm rotation', pi};

% Tracking markers
forearm.markers = {'styl_ulnae' styl_uln
		  'styl_radii' styl_rad
		  'radius' radius};

%-------------------
% palm
%-------------------

palm.name = 'palm';
palm_center = wjc; 
palm.localframe = cat(1, cat(2, e_x, e_y, e_z, wjc),...
		       [0 0 0 1]);
palm.dof = {[1 2], []}; %
palm.states = {'wrist flex', pi/2
	        'wrist adduction', pi/6};
palm.markers = { 'mcp2', mcp2
		 'mcp3', mcp3
		 'mcp4', mcp4
		 'mcp5', mcp5};


%-------------------
% Thumb base
%-------------------


thumbb.name = 'thumb base';
thumbb.center = styl_rad; 
e_z = mcp1 - styl_rad;
e_z = e_z / norm(e_z);
e_x = e_ml;
e_x = e_x - (e_x'*e_z)*e_z;
e_x = e_z / norm(e_x);
e_y = cross(e_z, e_x);

thumbb.localframe = cat(1, cat(2, e_x, e_y, e_z, thumbb.center),...
		       [0 0 0 1]);
thumbb.dof = {[1 2], []}; %
thumbb.states = {'thumb flex', pi/3
	        'thumb adduction', pi/3};
thumbb.markers = {'mcp1', mcp1};


%-------------------
% Thumb dp
%-------------------

thumbdp.name = 'thumb dp';
e_z = tip1 - pip1;
e_z = e_z / norm(e_z);
e_x = e_ml;
e_x = e_x - (e_x'*e_z)*e_z;
e_x = e_x / norm(e_x);
e_y = cross(e_z, e_x);
thumbdp.center = pip1 + 0.005*e_y; % 5mm below pip marker 
thumbdp.CoM = (0.75*tip1 + 0.25*pip1) + 0.005*e_y; % 10mm below pip marker 
thumbdp.radius = 0.01;

%%thumbdp.CoM = plate1; %% The point of contact with the plate

thumbdp.localframe = cat(1, cat(2, e_x, e_y, e_z, thumbdp.center),...
		       [0 0 0 1]);
thumbdp.dof = {[1], []}; %
thumbdp.states = {'pip1 flex', pi/2};
thumbdp.markers = {'tip1', tip1};


[tws2, p02, gcnames2, jc2, segmnames2, CoM2, rad2] = ...
    build_finger('2', mcp2, dip2, pip2, tip2, mcp3, plate2);
[tws3, p03, gcnames3, jc3, segmnames3, CoM3, rad3] = ...
    build_finger('3', mcp3, dip3, pip3, tip3, mcp4, plate3);
[tws4, p04, gcnames4, jc4, segmnames4, CoM4, rad4] = ...
    build_finger('4', mcp4, dip4, pip4, tip4, mcp5, plate4);
[tws5, p05, gcnames5, jc5, segmnames5, CoM5, rad5] = ...
    build_finger('5', mcp5, dip5, pip5, tip5, mcp4, plate5);


%----------------------------------------------------------------
% Define the complete model
%----------------------------------------------------------------

[tws, p0, gcnames, jc, segmnames, CoM, radius] = build_model(forearm);
[tws_p, p0_p, gcnames_p, jc_p, segmnames_p, CoM_p, radius_p] = ...
    build_model(palm);
[tws1, p01, gcnames1, jc1, segmnames1, CoM1, rad1] = ...
    build_model(thumbb, thumbdp);

tws_p{2} = tws1;
tws_p{3} = tws2;
tws_p{4} = tws3;
tws_p{5} = tws4;
tws_p{6} = tws5;

tws{2} = tws_p;

p0_p{2} = p01;
p0_p{3} = p02;
p0_p{4} = p03;
p0_p{5} = p04;
p0_p{6} = p05;

p0{2} = p0_p;

jc_p{2} = jc1;
jc_p{3} = jc2;
jc_p{4} = jc3;
jc_p{5} = jc4;
jc_p{6} = jc5;

jc{2} = jc_p;

gcnames_p = cat(1, gcnames_p, ...
		gcnames1, ...
		gcnames2, ...
		gcnames3, ...
		gcnames4, ...
		gcnames5); 
gcnames = cat(1, gcnames, gcnames_p);

segmnames_p = cat(1, segmnames_p, ...
		  segmnames1, ...
		  segmnames2, ...
		  segmnames3, ...
		  segmnames4, ...
		  segmnames5); 
segmnames = cat(1, segmnames, segmnames_p);

CoM_p{2} = CoM1;
CoM_p{3} = CoM2
CoM_p{4} = CoM3;
CoM_p{5} = CoM4;
CoM_p{6} = CoM5;
CoM{2} = CoM_p;

radius_p{2} = rad1;
radius_p{3} = rad2;
radius_p{4} = rad3;
radius_p{5} = rad4;
radius_p{6} = rad5;
radius{2} = radius_p;

tm.twists = tws;
tm.p0 = p0;
tm.jcs = jc;
tm.gcnames = gcnames;
tm.segm_names = segmnames;
tm.CoM = CoM;
tm.radius = radius;

function m = myextractmeanmarkers(rd, mname)
% Will look in struct rd for marker (or landmark) of name
% mname. Returns the average position in a 3x1 column vector, or
% NaNs if not found.

m = nan(3,1);

if isstruct(rd)
  
  if isfield(rd,mname)
    md = getfield(rd,mname);
    keyboard
    m = (mean(md{1}(:,1:3),1))';
  end

else
  mm = extractmarkers(rd,mname);
  m = mm(1,1:3)';
end

endfunction

function [tws_f, p0_f, gcnames_f, jc_f, segmnames_f, CoM_f, radius_f] = \
      build_finger(name, mcp2, dip2, pip2, tip2, mcp3, plate);

%-------------------
% mc
%-------------------

mc2.name = ['metacarpal', name];
e_x = mcp3 - mcp2; % vector between neighbouring mcp markers. Pointing lateral 
e_x = e_x / norm(e_x);
e_z = dip2 - pip2;
e_z = e_z - (e_z'*e_x)*e_x;
e_z = e_z / norm(e_z);
e_y = cross(e_z, e_x);
mc2.center = mcp2 + 0.005*e_y; % 5mm below mcp marker 

mc2.localframe = cat(1, cat(2, e_x, e_y, e_z, mc2.center),...
		       [0 0 0 1]);
mc2.dof = {[1 2], []}; %
mc2.states = {['mcp', name, ' f/e'], pi/2
	      ['mcp', name, ' a/a'], pi/4};
mc2.markers = {['pip', name], pip2};


%-------------------
% pp2
%-------------------

pp2.name = ['proxphal', name];
e_z = dip2 - pip2;
e_z = e_z - (e_z'*e_x)*e_x;
e_z = e_z / norm(e_z);
e_y = cross(e_z, e_x);
pp2.center = pip2 + 0.004*e_y; % 4mm below pip marker 

pp2.localframe = cat(1, cat(2, e_x, e_y, e_z, pp2.center),...
		       [0 0 0 1]);
pp2.dof = {[1], []}; %
pp2.states = {['pip', name, ' f/e'], pi/2};
pp2.markers = {['dip', name], dip2};



%-------------------
% dp2
%-------------------

dp2.name = ['distphal', name];
e_z = tip2 - dip2;
e_z = e_z - (e_z'*e_x)*e_x;
e_z = e_z / norm(e_z);
e_y = cross(e_z, e_x);
dp2.center = dip2 + 0.003*e_y; % 3mm below dip marker 
dp2.CoM = (0.75*tip2 + 0.25*dip2) + 0.005*e_y; % 10mm below pip marker 
dp2.radius = 0.01;
%%dp2.CoM = tip2 + 0.006*e_y; % Used for endpoing 6mm below tip marker
dp2.CoM = plate; % The point of contact with the surface
dp2.localframe = cat(1, cat(2, e_x, e_y, e_z, dp2.center),...
		       [0 0 0 1]);
dp2.dof = {[1], []}; %
dp2.states = {['dip', name, ' f/e'], pi/2};
dp2.markers = {['tip', name], tip2};

[tws_f, p0_f, gcnames_f, jc_f, segmnames_f, CoM_f, radius_f] = ...
    build_model(mc2, pp2, dp2);


endfunction