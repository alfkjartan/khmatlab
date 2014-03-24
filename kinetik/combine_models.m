function m = combine_models(mroot, m1, m2)
%%  m = combine_models(mroot, m1, m2)
%% Will make m1 and m2 branches of model mroot. If mroot is empty, then a simple model
%% with no degrees of freedom and no mass will be generated to hold the two branches.

%% Kjartan Halvorsen
%% 2013-07-16

if isempty(mroot)
  mroot.name = 'mroot';
  mroot.localframe = [eye(3), zeros(3,1)
		      zeros(1,3), 1];
  mroot.dof = {[],[]};
  mroot.states = {};
  mroot.CoM = zeros(3,1);
  mroot.g0 = [eye(3) mroot.CoM
	      zeros(1,3), 1];
  mroot.mass = 0;
  mroot.moment_of_inertia = 0*eye(3);
  mroot.generalized_inertia = [0*eye(3), zeros(3,3)
			       zeros(3,3), mroot.moment_of_inertia];

  [tws, p0, gcnames, jc, segmnames, CoM, radius, mass, g0, inertia] = build_model(mroot);
  mroot.twists = tws;
  mroot.p0 = p0;
  mroot.gcnames = gcnames;
  mroot.jcs = jc;
  mroot.segm_names = segmnames;
  mroot.CoM = CoM;
  mroot.mass = mass;
  mroot.g0 = g0;
  mroot.inertia = inertia;
end

m = mroot;
m.twists{2} = m1.twists;
m.twists{3} = m2.twists;

m.p0{2} = m1.p0;
m.p0{3} = m2.p0;

m.jcs{2} = m1.jcs;
m.jcs{3} = m2.jcs;

m.gcnames{2} = m1.gcnames;
m.gcnames{3} = m2.gcnames;

m.segm_names{2} = m1.gcnames;
m.segm_names{3} = m2.gcnames;

m.CoM{2} = m1.CoM;
m.CoM{3} = m2.CoM;

m.g0{2} = m1.g0;
m.g0{3} = m2.g0;

m.inertia{2} = m1.inertia;
m.inertia{3} = m2.inertia;

m.mass{2} = m1.mass;
m.mass{3} = m2.mass;

