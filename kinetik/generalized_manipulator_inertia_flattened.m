function [M, Mtrue] = generalized_manipulator_inertia_flattened(Mb, ...
                                                  tws, g0, theta)
%  M = generalized_manipulator_inertia_flattened(Mb, ...
%                                                  tws, g0, theta)
% Returns the generalized inertia for the manipulator described by the model struct km.
% See eq (4.29) in Murray, Li, Sastry 
%
%% Input
%    Mb         ->  The local inertia matrices (6 x 6 x ndofs)
%    tws        ->  The twists (6 x ndofs)
%    g0         ->  the transformations from local to spatial frame
%                   of the base for theta=0 (4 x 4 x ndofs)
%    theta      ->  The joint angles, (ndofs x 1)
%% Output
%    M          <-  The generalized manipulator inertia matrix
%                   (ndofs x ndofs ) i.e. the inertia matrix in joint space.

%% Kjartan Halvorsen
% 2017-03-20

if nargin == 0
   [M, Mtrue] = do_unit_test();
else

    ndofs = size(theta, 1);
    
    Mlinks = link_inertia_flattened(Mb, g0)
    
    M = zeros(ndofs, ndofs);

    for i=1:ndofs
        for j=1:ndofs
            for l = max(i,j):ndofs
                Ali = link_adjoint_trf(tws, theta, l, i);
                Alj = link_adjoint_trf(tws, theta, l, j);
                M(i,j) = M(i,j) + ...
                         tws(:,i)'*Ali'*Mlinks(:,:,l)*Alj*tws(:,j);
            end
        end
    end
end

end


function [M, Mtrue] = do_unit_test()

%% Construct model of the scara robot
l1 = 1;
l2 = 2;
m1 = 1;
m2 = 2;
m3 = 0.6;
m4 = 0.3;

sm = scara_robot_model(l1, l2, m1, m2, m3, m4);

states = [pi/4; pi/5; pi/6; pi/7];

Mtrue0 = scara_robot_inertia(sm, states);
Mtrue = generalized_manipulator_inertia(sm, states);

[tws, g0, Mb] = flatten_km(sm);
M = generalized_manipulator_inertia_flattened(Mb, tws, g0, states);

tol = 1e-12;

if norm(Mtrue - Mtrue0) > tol
   disp('Test 1 failed')
   disp('Expected'), Mtrue
   disp('Found'), M
else
   disp('Test 1 OK')
end

if norm(M - Mtrue) > tol
   disp('Test 2 failed')
   disp('Expected'), Mtrue
   disp('Found'), M
else
   disp('Test 2 OK')
end

if norm(M - Mtrue0) > tol
   disp('Test 3 failed')
   disp('Expected'), Mtrue0
   disp('Found'), M
else
   disp('Test 3 OK')
end

end