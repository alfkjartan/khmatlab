function [tm, endpointstr] = build_simple_model(refdata)
%  gm = build_simple_model(refdata)
%
% Returns a kinematic model for the hand tracking data
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



[tws, p0, gcnames, jc, segmnames, CoM, radius] = build_model(forearm);
[tws_p, p0_p, gcnames_p, jc_p, segmnames_p, CoM_p, radius_p] = ...
    build_model(palm);
tws{2} = tws_p;
p0{2} = p0_p;
jc{2} = jc_p;
gcnames = cat(1, gcnames, gcnames_p);
segmnames = cat(1, segmnames, segmnames_p);
CoM{2} = CoM_p;
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