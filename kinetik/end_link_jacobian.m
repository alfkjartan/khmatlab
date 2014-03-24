function JJ = end_link_jacobian(km, states)
%%  JJ = end_link_jacobian(km, states)
%% Returns the jacobian of the end link, defined by 
%%    JJ*\dot{states} = \dot{x},
%% where \dot{x} is the velocity (instantanous twist) of the end link in body coordinates.
%% If the kinematic model contains branches, JJ will have a block column structure, with
%% the jacobians for the different end points stacked.
%%
%% Input
%%    km         ->  kinematic model struct
%%    states     ->  joint states (nsts x nfrs)
%% Output
%%    JJ         <-  Jacobian (JJ1 ; JJ2 ; ... ; JJm) for m branches (m endpoints),
%%                   where JJi is (6 x nsts x nfrs)


%% Kjartan Halvorsen
%% 2013-06-21

if (nargin == 0)
   do_unit_test();
   return
end

%% Get twists and corresponding states for each branch
[branchtws, branchg0, branchstates, brinds] = get_branches(km, states); 

[nsts, nfrs] = size(states);

JJi = body_manipulator_jacobian(km.twists, km.g0, states(:,1));
nbrs = size(JJi,3);

JJ = zeros(6*nbrs, nsts, nfrs);

for  i=1:nfrs
  JJi = body_manipulator_jacobian(km.twists, km.g0, states(:,i));
  for br = 1:nbrs
    JJ(6*(br-1)+1:6*br, :, i) = JJi(:,:,br);
  end
end


function do_unit_test()
	 km1 = planar_link_model();
	 km2 = planar_link_model();

	 km.twists = {{}, km1.twists, km2.twists};
	 km.g0 = {eye(4,4), km1.g0, km2.g0};

	 nsts = 6;
	 nfrs = 3;
	 states = randn(nsts, nfrs);

	 JJ = end_link_jacobian(km, states);
	 size(JJ)
	 JJ1 = zeros(6, nsts, nfrs);
	 JJ2 = zeros(6, nsts, nfrs);
	 JJbmj = zeros(12, nsts, nfrs);
	 for i = 1:nfrs
	     JJ1i = link_jacobians(km1.twists, km1.g0, states([1 2 3],i));
	     JJ2i = link_jacobians(km2.twists, km2.g0, states([4 5 6],i));
	     JJ1(:,[1 2 3],i) = JJ1i(:,:,end); 
	     JJ2(:,[4 5 6],i) = JJ2i(:,:,end); 
	     JJbmji = body_manipulator_jacobian(km.twists, km.g0, states(:,i));
	     JJbmj(1:6, :, i) = JJbmji(:,:,1);
	     JJbmj(7:12, :, i) = JJbmji(:,:,2);
	 end

	 assert(JJ(1:6, :, :), JJ1, 1e-12);
	 disp('Test 1 OK')
	 assert(JJ(7:12, :, :), JJ2, 1e-12);
	 disp('Test 2 OK')

	 keyboard
	 assert(JJ, JJbmj, 1e-12);
	 disp('Test 3 OK')

