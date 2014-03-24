% Script that loads clubmarker and clublandmark locations (exported from
% Visual3D), and computes CoM location and inertia matrix for the combined
% club segment consisting of epoxygrip, clubshaft, hosel and clubhead. 
%
% Fredrik Tinmark 2013-03-27

% Segment masses
mass_ch = 0.235;
mass_cs = 0.119;
mass_eg = 0.104;
mass_h = 0.067;

club_mass = mass_ch + mass_cs + mass_eg + mass_h;

% Moment of inertias for the individual club segments about their principle
% axes (from Visual3D)
I_principal_ch = [6.20045e-005 0 0; 0 7.80155e-005 0; 0 0 2.43018e-005];
% I_principal_ch = [7.08197e-005 0 0; 0 7.89414e-005 0; 0 0 1.64126e-005];
I_principal_cs = [0.0065626 0 0; 0 0.0065626 0; 0 0 4.84968e-006];
I_principal_eg = [0.000812585 0 0; 0 0.000816024 0; 0 0 5.2196e-006];
I_principal_h = [2.93544e-005 0 0; 0 2.93544e-005 0; 0 0 1.53767e-006];

refdata = load('C:\Users\fredrikt\Documents\MATLAB\PlianceValidation\20120223\clubref.mat');

clubhead_CoM = mmyextractmeanmarkers(refdata, 'CoM_ClubHead');
clubshaft_CoM = mmyextractmeanmarkers(refdata, 'CoM_ClubShaft');
epoxygrip_CoM = mmyextractmeanmarkers(refdata, 'CoM_EpoxyGrip');
hosel_CoM = mmyextractmeanmarkers(refdata, 'CoM_Hosel');

clubheel_bottom = mmyextractmeanmarkers(refdata, 'CHB');
clubheel_top = mmyextractmeanmarkers(refdata, 'CHT');
clubtoe_bottom = mmyextractmeanmarkers(refdata, 'CTB');
clubtoe_top = mmyextractmeanmarkers(refdata, 'CTT');

club_top = mmyextractmeanmarkers(refdata, 'ClubTop');
epoxygrip_bottom = mmyextractmeanmarkers(refdata, 'EpoxyGripBottom');
hosel_distal = mmyextractmeanmarkers(refdata, 'HoselDistal');
shaft_bottom = mmyextractmeanmarkers(refdata, 'ShaftBottom');

% CoM coordinates for the combined club in the global frame
club_CoM = (mass_ch*clubhead_CoM + mass_cs*clubshaft_CoM + ...
    mass_eg*epoxygrip_CoM + mass_h*hosel_CoM) / club_mass;

% Global frame
LAB = eye(3);

% Local frame for the combined club segment (parallel to the axes of the
% epoxy grip)
e_IS = club_top - epoxygrip_bottom;
e_IS = e_IS / norm(e_IS);
e_PA = clubtoe_bottom - clubheel_bottom;
e_PA = e_PA - (e_PA'*e_IS)*e_IS; 
e_PA = e_PA / norm(e_PA);
e_LR = cross(e_PA, e_IS);

% Rotation matrix that relates the orientation of the global frame to the
% combined club frame
R = [e_LR, e_PA, e_IS]'*LAB;

% CoM coordinates for the individual club segmants in the combined club
% frame
com_ch = R*(clubhead_CoM - club_CoM);
com_cs = R*(clubshaft_CoM - club_CoM);
com_h = R*(hosel_CoM - club_CoM);
com_eg = R*(epoxygrip_CoM - club_CoM);

% Local frame for the clubhead segment
e_z = (0.5*clubheel_top + 0.5*clubheel_bottom) - (0.5*clubtoe_top + 0.5*clubtoe_bottom);
e_z = e_z / norm(e_z);
e_x = clubtoe_top - clubtoe_bottom;
e_x = e_x - (e_x'*e_z)*e_z;
e_x = e_x / norm(e_x);
e_y = cross(e_z, e_x);

% Rotation matrix that relates the orientation of the local frame for the
% clubhead segment to the combined club frame
R = [e_LR, e_PA, e_IS]'*[e_x, e_y, e_z];

% Inertia for the clubhead segment with respect to the combined club frame
% with origin at the club CoM
[I_chl, I_chr] = golf_club_inertia_tensor(R, I_principal_ch, mass_ch, com_ch);

% Local frame for the clubshaft segment
e_z = club_top - shaft_bottom;
e_z = e_z / norm(e_z);
e_y = clubtoe_bottom - clubheel_bottom;
e_y = e_y - (e_y'*e_z)*e_z; 
e_y = e_y / norm(e_y);
e_x = cross(e_y, e_z);

% Rotation matrix that relates the orientation of the local frame for the
% clubshaft segment to the combined club frame
R = [e_LR, e_PA, e_IS]'*[e_x, e_y, e_z];

% Inertia for the clubshaft segment with respect to the combined club frame
[I_csl, I_csr] = golf_club_inertia_tensor(R, I_principal_cs, mass_cs, com_cs);

% Local frame for the hosel segment
e_z = shaft_bottom - hosel_distal;
e_z = e_z / norm(e_z);
e_y = clubtoe_bottom - clubheel_bottom;
e_y = e_y - (e_y'*e_z)*e_z; 
e_y = e_y / norm(e_y);
e_x = cross(e_y, e_z);

% Rotation matrix that relates the orientation of the local frame for the
% hosel segment to the combined club frame
R = [e_LR, e_PA, e_IS]'*[e_x, e_y, e_z];

% Inertia for the hosel segment with respect to the combined club frame
[I_hl, I_hr] = golf_club_inertia_tensor(R, I_principal_h, mass_h, com_h);

% Local frame for the epoxygrip segment
e_z = e_IS;
e_y = e_PA;
e_x = e_LR;

R = eye(3);

% Inertia for the epoxygrip segment with respect to the combined club frame
[I_egl, I_egr] = golf_club_inertia_tensor(R, I_principal_eg, mass_eg, com_eg);

% Inertia for the combined club segment
I_l = I_chl + I_csl + I_egl + I_hl;
I_r = I_chr + I_csr + I_egr + I_hr;

% The principal axes and principal moments of inertia for the combined club segment
I_cc = sparse(I_l + I_r);
[V,D] = eigs(I_cc);
W = V'*[-1 0 0; 0 1 0; 0 0 1];

% The coordinates of the CoM in a frame parallel to principal axes
% with origin at ClubTop
CoM_local = W*[e_x, e_y, e_z]'*(club_CoM - club_top);
