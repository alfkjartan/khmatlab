function [mass, CoM, local_inertia, inertia] = get_club_model(grip_top, toe_top_grove, toe_bottom_grove, heel_bottom_grove)
%%  [mass, CoM, local_inertia, inertia] = get_club_model(grip_top, toe_top_grove, toe_bottom_grove, heel_bottom_grove)
%%  Returns a model of the club (inertial properties). Masses are hardcoded below
%%
%%  Input arguments
%%     grip_top       ->    position of point on top of grip
%%     toe_top_grove  ->    position of point on toe end of top grove
%%     toe_bottom_grove  ->    position of point on toe end of bottom grove
%%     heel_bottom_grove  ->    position of point on heel end of bottom grove
%%  Output
%%     mass       <-   mass of the complete club
%%     CoM        <-   Center of mass of the complete club
%%     local_inertia  <-  moment of inertia tensor wrt CoM in local coordinate system
%%                        of the club (z points down the shaft, y points along the bottom grove,
%%                        x is perpendicular to z and y)
%%     inertia        <-   moment of inertia tensor wrt CoM in global coordinate system.

%% Kjartan Halvorsen 
%% 2013-10-09


%% The club is assumed to be made up of three parts: grip, shaft and head. 
%% The moment of inertia is determined for each of these parts, and combined
%% to form a hand+club segment. The steps are:
%%   1) Determine inertia matrix in local coordinate system of part.
%%   2) Transfrom to local coordinate system of hand
%%   3) Add to inertia of hand using the paralell axis theorem.

%% Hardcoded parameters
grip.mass = 50*1e-3; % 50 g
grip.length = 0.27; % 27cm typical length of the grip
grip.r2 = 0.01; % Radius of outer contour of the grip
grip.r1 = 0.006; % 6mm radius of shaft 
shaft.mass = 120*1e-3; % 120 g
shaft.r = 0.004; % average radius of  shaft
clubhead.mass = 300*1e-3; % 300 g

e_shaft = heel_bottom_grove - grip_top; % Points down the shaft
e_shaft = e_shaft / norm(e_shaft);
e_club_pa = toe_bottom_grove - heel_bottom_grove;
e_club_pa = e_club_pa - (e_club_pa'*e_shaft)*e_shaft;
e_club_pa = e_club_pa / norm(e_club_pa);
e_club_x = cross(e_club_pa, e_shaft);
R_club = [e_club_x e_club_pa e_shaft];

club_head_center = 0.5*heel_bottom_grove + 0.5*toe_top_grove; 

grip.CoM = grip_top + 0.4*grip.length*e_shaft; % CoM 12cm down the shaft. Grip typically 27cm long.
grip.local_inertia = grip.mass * diag( [1/12*(3*(grip.r1^2 + grip.r2^2) + grip.length*2)
					1/12*(3*(grip.r1^2 + grip.r2^2) + grip.length*2)
					1/2 * (grip.r1^2 + grip.r2^2)]);
grip.inertia = R_club * grip.local_inertia * R_club';

shaft.length = norm(grip_top - heel_bottom_grove); % length of the shaft
shaft.CoM = grip_top + 0.35*shaft.length*e_shaft; 
shaft.local_inertia = shaft.mass * diag( [1/12*(6*shaft.r^2 + shaft.length*2)
					  1/12*(6*shaft.r^2 + shaft.length*2)
					  shaft.r^2]);
shaft.inertia = R_club * shaft.local_inertia * R_club';

%% Clubhead as point mass
clubhead.CoM = (0.5* toe_bottom_grove + 0.5*heel_bottom_grove);
clubhead.inertia = zeros(3,3);

%% Combine
[inertia, CoM, mass] = combine_inertia(grip.inertia, grip.CoM, grip.mass, shaft.inertia, shaft.CoM, shaft.mass);
[inertia, CoM, mass] = combine_inertia(inertia, CoM, mass, clubhead.inertia, clubhead.CoM, clubhead.mass);

local_inertia = R_club' * inertia * R_club;

