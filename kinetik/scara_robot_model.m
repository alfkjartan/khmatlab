function sm = scara_robot_model(l1, l2, m1, m2, m3, m4, pos0, dirarm)
% Returns a model of the scara robot with link lengths l1, l2 and masses m1, m2, m3, m4

%% Kjartan Halvorsen
% 2013-06-05

if nargin == 0
   % Defaults
   l1 = 1;
   l2 = 2;
   m1 = 1;
   m2 = 2;
   m3 = 0.6;
   m4 = 0.3;
   pos0 = zeros(3,1);
   dirarm = [0;1;0]; %% Direction of the arm
end

if nargin < 8
   dirarm = [0;1;0];
end

if nargin < 7
   pos0 = zeros(3,1);
end


Z3 = zeros(3,3);
I3 = eye(3);

link1.name = 'link1';
link1.localframe = [I3, pos0
		    zeros(1,3), 1];
link1.dof = {[3],[]};
link1.states = {'th1', pi};
link1.CoM = pos0 + l1/2*dirarm;
ey = dirarm;
ez = [0;0;1];
ex = cross(ey,ez);
link1.g0 = [ex ey ez link1.CoM
	   zeros(1,3), 1];
link1.mass = m1;
link1.moment_of_inertia = m1 * diag( [l1^2/12; 0; l1^2/12] );
link1.generalized_inertia = [m1*I3, Z3
			     Z3, link1.moment_of_inertia];

link2.name = 'link2';
link2.localframe = [I3, pos0+l1*dirarm
		    zeros(1,3), 1];
link2.dof = {[3],[]};
link2.states = {'th2', pi};
link2.CoM = pos0 + (l1+l2/2)*dirarm;
link2.g0 = [ex ey ez, link2.CoM
	   zeros(1,3), 1];
link2.mass = m2;
link2.moment_of_inertia = m2 * diag( [l2^2/12; 0; l2^2/12] );
link2.generalized_inertia = [m2*I3, Z3
			     Z3, link2.moment_of_inertia];

link3.name = 'link3';
link3.localframe = [I3 pos0 + (l1+l2)*dirarm
		    zeros(1,3), 1];
link3.dof = {[3],[]};
link3.states = {'th3', pi};
link3.CoM = pos0 + (l1+l2)*dirarm;
link3.g0 = [ex ey ez link3.CoM
	   zeros(1,3), 1];
link3.mass = m3;
link3.moment_of_inertia = m3 * diag( [l1^2/12; l1^2/12; l2^2/24] );
link3.generalized_inertia = [m3*I3, Z3
			     Z3, link3.moment_of_inertia];

link4.name = 'endlink';
link4.localframe = [I3, pos0 + (l1+l2)*dirarm
		    zeros(1,3), 1];
link4.dof = {[],[3]};
link4.states = {'th4', 0.5};
link4.CoM = pos0 + (l1+l2)*dirarm;
link4.g0 = [I3, link4.CoM
	   zeros(1,3), 1];
link4.mass = m4;
link4.moment_of_inertia = m4 * diag( [l1^2/12; l1^2/12; l2^2/24] );
link4.generalized_inertia = [m4*I3, Z3
			     Z3, link4.moment_of_inertia];

[tws, p0, gcnames, jc, segmnames, CoM, radius, mass, g0, inertia] = build_model(link1, link2, link3, link4);

sm.twists = tws;
sm.p0 = p0;
sm.jcs = jc;
sm.gcnames = gcnames;
sm.segm_names = segmnames;
sm.CoM = CoM;
sm.g0 = g0;
sm.inertia = inertia;
sm.mass = mass;
sm.l1 = l1;
sm.l2 = l2;

