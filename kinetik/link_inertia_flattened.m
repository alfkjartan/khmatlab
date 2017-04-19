function Ms = link_inertia_flattened(Mb, gsli)
%   Ms = link_inertia_flattened(Mb, gsli)
% Returns a list of link inertias for the flattened kinematic chain
% (by flattened means the link is defined with as many links as
% degrees of freedom, some links have zero length and zero mass). 
%
% Input
%   Mb    ->   Inertia matrices for each link in local coordinate
%              system (6 x 6 x ndofs)
%   gsli  ->   Transformations from each link's local coordinate
%              system to the spatial coordinate system of the base.
%
% Output
%   Ms     <-  Inertia matrices reflected into the coordinate sytem
%              of the base

% Kjartan Halvorsen


if nargin == 0
   do_unit_test();
else
   tw0 = zeros(6,1);
   Ms = zeros(size(Mb));
   for i=1:size(Mb, 3)
       [twslask, Ad_inv] = adjoint_trf_inv(tw0, gsli(:,:,i));
       Ms(:,:,i) = Ad_inv'*Mb(:,:,i)*Ad_inv;
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
m = [m1;m2;m3;m4];

states = [pi/4; pi/5; pi/6; pi/7];

sm = scara_robot_model(l1, l2, m1, m2, m3, m4);
[M, M1, M2, M3, M4] = scara_robot_inertia(sm, states)

MM = cat(3, M1, M2, M3, M4);

Ms0 = zeros(6,6,4);
CoMs = cat(2, [0;0.5*l1;0],...
           [0; l1+0.5*l2; 0],...
           [0; l1 + l2; 0],...
           [0; l1+l2; 0]);

Adjs0 = zeros(6,6,4);
Adjs = zeros(6,6,4);
for i=1:4
    p = CoMs(:,i);
    phat = hat(p);
    gi = eye(4);
    gi(1:3, 4) = CoMs(:,i);
    
    [slsk, Adjs(:,:,i)] = adjoint_trf_inv(zeros(6,1), gi);
    % From page 177 , corrected
    Adjs0(:,:,i) = eye(6);
    Adjs0(1:3,4:6,i) = -phat;
    Ms0(1:3,1:3,i) = eye(3)*m(i);
    Ms0(4:6,4:6,i) = MM(4:6, 4:6, i) - m(i)*phat*phat;
    Ms0(1:3, 4:6, i) = -m(i)*phat;
    Ms0(4:6, 1:3, i) = m(i)*phat;
end

[tws, g0, Mb] = flatten_km(sm);

Ms = link_inertia_flattened(Mb, g0);

tol = 1e-12;

for i=1:4
    if norm(Adjs(:,:,i) - Adjs0(:,:,i)) > tol
        disp(sprintf('Test %d a) failed', i))
        disp('Expected'), Adjs0(:,:,i)
        disp('Found'), Adjs(:,:,i)
    else
        disp(sprintf('Test %d a) OK', i))
    end
    if norm(Ms(:,:,i) - Ms0(:,:,i)) > tol
        disp(sprintf('Test %d b) failed', i))
        disp('Expected'), Ms0(:,:,i)
        disp('Found'), Ms(:,:,i)
    else
        disp(sprintf('Test %d b) OK', i))
    end
end

end