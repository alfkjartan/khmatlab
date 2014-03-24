function gm = build_golf_test_model()
%% Returns a test model. Basically one link branching into two planar linkages
%% with a common frame for the last link.

%% Kjartan Halvorsen
%% 2013-06-14

lengths = [1;1];
masses = [1;1];
posL = [0;0.5;0];
[armL, linksL] = planar_link_model(lengths, masses, posL);
posR = [0;-0.5;0];
[armR, linksR] = planar_link_model(lengths, masses, posR);


mb = 2; % mass of base link
lb = 1; % length of base link

I3 = eye(3);
Z3 = zeros(3,3);

baslink.name = 'base';
baslink.localframe = [I3, zeros(3,1)
		       zeros(1,3), 1];
baslink.dof = {[3],[]};
baslink.states = {sprintf('th%d', i), pi};
baslink.CoM = [0; 0; 0];
baslink.g0 = [I3, baslink.CoM
	       zeros(1,3), 1];
baslink.mass = mb;
baslink.moment_of_inertia = mb * diag( [lb^2/12; 0; lb^2/12] );
baslink.generalized_inertia = [mb*I3, Z3
				  Z3, baslink.moment_of_inertia];


[tws, p0, gcnames, jc, segmnames, CoM, radius, mass, g0, inertia] = build_model(baslink);
[twsL, p0L, gcnamesL, jcL, segmnamesL, CoML, radiusL, massL, g0L, inertiaL] = build_model(linksL{:});
[twsR, p0R, gcnamesR, jcR, segmnamesR, CoMR, radiusR, massR, g0R, inertiaR] = build_model(linksR{:});

tws{2} = twsL;
tws{3} = twsR;

p0{2} = p0L;
p0{3} = p0R;

gcnames{2} = p0L;
gcnames{3} = p0R;


