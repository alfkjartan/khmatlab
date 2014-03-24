function Jsb=body_manipulator_jacobian(tws, gsli, th, twsprox, gprox)
%% Jsb=body_manipulator_jacobian(tws, gsli, th, twsprox, gprox)
%% Returns the body manipulator jacobian for the kinematic linkage
%% Input
%%   tws    ->  nested array of twists
%%   gsli   ->  nested array of transformations taking points in the local coordinate system of
%%              each link to the static frame for th=0
%%   th     ->  generalized coordinates
%%   twsprox -> twists from root up to current segment (6 x nprox) matrix
%%   gprox   -> proximal rigid transformations (4 x 4 x nprox) matrix, where the 
%%              ith transformation is
%%                 g_i = expr(twsprox(1), thprox(1))* ... * expr(twsprox(i), thprox(i))
%% Output
%%   Jsb    <-  The jacobians, a (6 x nsts x n_endpoints) matrix

%% Kjartan Halvorsen
%% 2013-06-04

if nargin == 0
%%   do_unit_test();
return
end

Jst = spatial_manipulator_jacobian(tws, th); % Returns a (6 x nsts x n_endpointsx) matrix

gst = forward_map(tws, gsli, th); %% returns a (4 x 4 x n_endpoints) matrix

n_endpoints = size(gst,3);

Jsb = zeros(size(Jst));
for i=1:n_endpoints
    Jsb(:,:,i) = adjoint_trf_inv( Jst(:,:,i), gst(:,:,i) ) ;
end
