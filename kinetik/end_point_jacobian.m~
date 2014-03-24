function JJ = end_point_jacobian(km, states)
%%  JJ = end_point_jacobian(km, states)
%% Returns the jacobian of the end point in spatial coordinates, defined by 
%%    JJ*\dot{states} = \dot{x},
%% where \dot{x} is the velocity in spatial coordinates of the end point.
%% If the kinematic model contains branches, JJ will have a block column structure, with
%% the jacobians for the different end points stacked.
%%
%% Input
%%    km         ->  kinematic model struct
%%    states     ->  joint states (nsts x nfrs)
%% Output
%%    JJ         <-  Jacobian (JJ1 ; JJ2 ; ... ; JJm) for m branches (m endpoints),
%%                   where JJi is (3 x nsts x nfrs)


%% Kjartan Halvorsen
%% 2013-07-08

if (nargin == 0)
   do_unit_test();
   return
end

[nsts, nfrs] = size(states);

JJs1 = spatial_manipulator_jacobian(km.twists, states(:,1)); % (6 x nsts x nbranches)
nbranches = size(JJs1,3);

JJ = zeros(3*nbranches, nsts, nfrs);

if isfield(km, 'object_frame')
  g0 = km.object_frame;
else
  g0 = km.g0;
end

%disp('### debug end_point_jacobian')
%disp('g0')
%g0

for i=1:nfrs
  JJsi = spatial_manipulator_jacobian(km.twists, states(:,i)); % (6 x nsts x nbranches)
  gsti = forward_map(km.twists, g0, states(:,i));

  for br = 1:nbranches
    %% The end point is defined as the origin of the tool frame (or body frame of the end 
    %% segment.
    pbr = gsti(:,4,br);
    for j = 1:nsts
	vj = hat(JJsi(:,j,br))*pbr;
	JJ((br-1)*3+1:br*3,j,i) = vj(1:3);
    end
  end
end


function do_unit_test()

	    l1 = 1;
	    l2 = 2;
	    m1 = 1;
	    m2 = 2;
	    m3 = 0.6;
	    m4 = 0.3;

	    sm = scara_robot_model(l1, l2, m1, m2, m3, m4);


	    delta = 1e-8;
	    th = randn(4,1);

	    JJ = end_point_jacobian(sm, th);

	    Jdif = zeros(3, 4);
	    gst = forward_map(sm.twists, sm.g0, th);
	    p0 = gst(1:3,4);
	    for st=1:4
	      thd = th;
	      thd(st) = thd(st) + delta;
	      gstd = forward_map(sm.twists, sm.g0, thd);
	      pd = gstd(1:3,4);
	      Jdif(:,st) = (pd-p0)/delta;
	    end

	    assert(JJ, Jdif, delta*10)
	    disp('Test 1 OK')
