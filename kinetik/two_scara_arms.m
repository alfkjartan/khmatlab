function sm = two_scara_arms()

  l1 = 1;
  l2 = 1;
  m1 = 1;
  m2 = 1;
  m3 = 1;
  m4 = 1;

  dirarm = [0;1;0];
  pos0 = -(l1+l2)*dirarm;
  sm1 = scara_robot_model(l1, l2, m1, m2, m3, m4, pos0, dirarm);
  dirarm = -[0;1;0];
  pos0 = -(l1+l2)*dirarm;
  sm2 = scara_robot_model(l1, l2, m1, m2, m3, m4, pos0, dirarm);

  link0.name = 'link0';
  link0.localframe = [eye(3), zeros(3,1)
		      zeros(1,3), 1];
  link0.dof = {[],[]};
  link0.states = {};
  link0.CoM = zeros(3,1);
  link0.g0 = [eye(3) link0.CoM
	      zeros(1,3), 1];
  link0.mass = 100;
  link0.moment_of_inertia = m1 * eye(3);
  link0.generalized_inertia = [m1*eye(3), zeros(3,3)
			       zeros(3,3), link0.moment_of_inertia];

  [tws, p0, gcnames, jc, segmnames, CoM, radius, mass, g0, inertia] = build_model(link0);

  tws{2} = sm1.twists;
  tws{3} = sm2.twists;
  sm.twists = tws;

  p0{2} = sm1.p0;
  p0{3} = sm2.p0;
  sm.p0 = p0;

  jc{2} = sm1.jcs;
  jc{3} = sm2.jcs;
  sm.jcs = jc;
  
  gcnames{2} = sm1.gcnames;
  gcnames{3} = sm2.gcnames;
  sm.gcnames = gcnames;

  segmnames{2} = sm1.gcnames;
  segmnames{3} = sm2.gcnames;
  sm.segm_names = segmnames;

  CoM{2} = sm1.CoM;
  CoM{3} = sm2.CoM;
  sm.CoM = CoM;

  g0{2} = sm1.g0;
  g0{3} = sm2.g0;
  sm.g0 = g0;

  inertia{2} = sm1.inertia;
  inertia{3} = sm2.inertia;
  sm.inertia = inertia;

  mass{2} = sm1.mass;
  mass{3} = sm2.mass;
  sm.mass = mass;
