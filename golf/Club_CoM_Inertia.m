mass_ch = 0.235;
mass_cs = 0.119;
mass_eg = 0.104;
mass_h = 0.067;

club_mass = mass_ch + mass_cs + mass_eg + mass_h;

refdata = load('C:\Users\fredrikt\Documents\MATLAB\PlianceValidation\20130411\clubref.mat');

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

LAB = eye(3);

% Local frame for the combined club segment (parallel to the axes of the
% epoxy grip)
e_IS = club_top - epoxygrip_bottom;
e_IS = e_IS / norm(e_IS);
e_PA = clubtoe_bottom - clubheel_bottom;
e_PA = e_PA - (e_PA'*e_IS)*e_IS; 
e_PA = e_PA / norm(e_PA);
e_LR = cross(e_PA, e_IS);

R = [e_LR, e_PA, e_IS]'*LAB;

% CoM coordinates for the combined club in the global frame
club_CoM = (mass_ch*clubhead_CoM + mass_cs*clubshaft_CoM + ...
    mass_eg*epoxygrip_CoM + mass_h*hosel_CoM) / club_mass;

v_ch = clubhead_CoM - club_CoM;
v_cs = clubshaft_CoM - club_CoM;
v_eg = epoxygrip_CoM - club_CoM;
v_h = hosel_CoM - club_CoM;

C = mass_ch*(v_ch*v_ch') + mass_cs*(v_cs*v_cs') + mass_eg*(v_eg*v_eg') + mass_h*(v_h*v_h');
I_g = eye(3)*trace(C) - C;

I_b = R*I_g*R';

v = club_CoM - club_top;

I_grepp_b = I_b + club_mass*(v'*v*eye(3) - v * v');

F_grepp = club_mass*(a - 9.81);

T_grepp = R'*I_grepp_b*R*alpha + cross(w ,R'*I_grepp_b*R)*w;