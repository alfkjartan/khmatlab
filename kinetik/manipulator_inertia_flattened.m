function C= manipulator_inertia_flattened(Mb, tws, g0, theta, ...
                                           thetadot, AA)
%  C = manipulator_inertia_flattened(Mb, tws, g0, theta, thetadot,
%                                    AA)
% Returns the Coriolis matrix for a flat kinematic chain.
% See eq (4.29) in Murray, Li, Sastry 
%
%% Input
%    Mb         ->  The local inertia matrices (6 x 6 x ndofs)
%    tws        ->  The twists (6 x ndofs)
%    g0         ->  the transformations from local to spatial frame
%                   of the base for theta=0 (4 x 4 x ndofs)
%    theta      ->  The joint angles, (ndofs x 1)
%    AA         ->  Optional, 6 x 6 x ndofs x ndofs array of
%                   adjoint transformations
%                   
%% Output
%    M          <-  The generalized manipulator inertia matrix
%                   (ndofs x ndofs ) i.e. the inertia matrix in joint space.

%% Kjartan Halvorsen
% 2017-03-28

if nargin == 0
   C = do_unit_test();
else

    ndofs = size(theta, 1);
    
    Mlinks = link_inertia_flattened(Mb, g0);
    
    if nargin < 6
        AA = [];
    end
    
    if isempty(AA)
        % Calculate all the adjoint transformations
        AA = zeros(6, 6, ndofs, ndofs);
        for i = 1:ndofs
            for j = 1:ndofs
                AA(:,:,i,j) = link_adjoint_trf(tws, theta, i, j);
            end
        end
    end
    
    
    % Calculate the inertia matrix
    M = zeros(ndofs, ndofs);
    for i=1:ndofs
        for j=1:ndofs
            M(i,j) = Mij(i,j, AA, tws, Mlinks);
        end
    end
    
end
end


function M = Mij(i,j,AA, tws, Mlinks)
% Computes the ij part of the inertia matrix. See MLS (4.29)
    
    ndofs = size(tws, 2);
    Mij = 0;

    for l = max(i,j):ndofs
        Mij = Mij + tws(:,i)' * AA(:,:,l,i)' * Mlinks(:,:,l) * AA(:,:,l,k) ...
             * tws(:,k);
    end
end
        
function C = do_unit_test()

%% Construct model of the scara robot
l1 = 1;
l2 = 2;
m1 = 1;
m2 = 2;
m3 = 0.6;
m4 = 0.3;

sm = scara_robot_model(l1, l2, m1, m2, m3, m4);

%theta = [pi/4; pi/5; pi/6; pi/7];
%thetadot = [pi/2; pi/1; pi/0.5; pi/0.25];
theta = randn(4,1);
thetadot = randn(4,1);

Mtrue = scara_robot_inertia(sm, theta, thetadot);


[tws, g0, Mb] = flatten_km(sm);
C = manipulator_coriolis_flattened(Mb, tws, g0, theta, thetadot);

tol = 1e-12;

if norm(C - Ctrue) > tol
   disp('Test 1 failed')
   disp('Expected'), Ctrue
   disp('Found'), C
else
   disp('Test 1 OK')
end

end