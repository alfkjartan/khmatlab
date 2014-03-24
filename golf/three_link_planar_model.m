function rm = three_link_planar_model(l1, l2, l3, m1, m2, m3)
%% Returns a model of a planar manipulator with link lenghts l1, l2, l3 and masses m1, m2, m3

%% Kjartan Halvorsen
%% 2013-06-12

if nargin == 0
   %% Defaults
   l1 = 1; l2 = 1; l3 = 1;
   m1 = 1; m2 = 1; m3 = 1;
end

Z3 = zeros(3,3);
I3 = eye(3);

link1.name = 'link1';
link1.localframe = [I3, zeros(3,1)
		    zeros(1,3), 1];
link1.dof = {[3],[]};
link1.states = {'th1', pi};
link1.CoM = [l1/2; 0; 0];
link1.g0 = [I3, link1.CoM
	   zeros(1,3), 1];
link1.mass = m1;
link1.moment_of_inertia = m1 * diag( [0; l1^2/12; l1^2/12] );
link1.generalized_inertia = [m1*I3, Z3
			     Z3, link1.moment_of_inertia];

link2.name = 'link2';
link2.localframe = [I3, [l1; 0;0]
		    zeros(1,3), 1];
link2.dof = {[3],[]};
link2.states = {'th2', pi};
link2.CoM = [l1+l2/2; 0; 0];
link2.g0 = [I3, link2.CoM
	   zeros(1,3), 1];
link2.mass = m2;
link2.moment_of_inertia = m2 * diag( [0; l2^2/12; l2^2/12] );
link2.generalized_inertia = [m2*I3, Z3
			     Z3, link2.moment_of_inertia];

link3.name = 'link3';
link3.localframe = [I3, [l1+l2;0;0]
		    zeros(1,3), 1];
link3.dof = {[3],[]};
link3.states = {'th3', pi};
link3.CoM = [l1+l2+l3/2; 0; 0];
link3.g0 = [I3, link3.CoM
	   zeros(1,3), 1];
link3.mass = m3;
link3.moment_of_inertia = m3 * diag( [0; l3^2/12; l3^2/12] );
link3.generalized_inertia = [m3*I3, Z3
			     Z3, link3.moment_of_inertia];

link3.p0 = [l1+l2+l3;0;0];
[tws, p0, gcnames, jc, segmnames, CoM, radius, mass, g0, inertia] = build_model(link1, link2, link3);


rm.twists = tws;
rm.p0 = p0;
rm.jcs = jc;
rm.gcnames = gcnames;
rm.segm_names = segmnames;
rm.CoM = CoM;
rm.g0 = g0;
rm.inertia = inertia;
rm.mass = mass;
rm.l1 = l1;
rm.l2 = l2;
rm.l3 = l3;

