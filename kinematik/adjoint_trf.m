function [twn, Ad_g] = adjoint_trf(tw, g)
%%  [twn, Ad_g] = adjoint_trf(tw, g)
%% Transforms the twist tw into twn by applying the adjoint transformation associated with g.
%% See eq 2.58 in Murray, Li, Sastry:
%%   Ad_g = [R  \hat{p}R ]
%%          [0      R    ],
%% where 
%%   g = [R  p]
%%       [0  1].

%% Kjartan Halvorsen
%% 2013-05-31

if nargin == 0
   do_unit_test();
else

  if size(tw,1) == 4
    twv = vee(tw);
  else
    twv = tw;
  end

  R = g(1:3, 1:3);
  p = g(1:3, 4);

  Ad_g = [ R  hat(p)*R
	   zeros(3,3) R];

  twn = Ad_g*twv;

  if size(tw, 1) == 4
    twn = hat(twn);
  end
end


function do_unit_test()

thr = 1e-12;

%% Check first that a twist is unchainged by the adjoint transformation corresponding to its
%% own exponential.

w1 = [0;0;1];
v1 = [1;0;0];
tw1 = [v1;w1];

g1 = expr(tw1, pi/4);

[twn1, Ad_g1] = adjoint_trf(tw1, g1);

if norm(tw1-twn1) > thr
   disp('Test 1 failed.')
   disp('Expected  '), tw1
   disp('Found  '), twn1
else
   disp('Test 1 OK.')
end

%% Check that inverse adjoint is same as adjoint of inverse

[twn1, Ad_g1_inv] = adjoint_trf_inv(tw1, g1);

g1inv = ginv(g1);
[twn1, Ad_g1inv] = adjoint_trf(tw1, g1inv);


if norm(Ad_g1_inv-Ad_g1inv) > thr
   disp('Test 2 failed.')
   disp('Expected  '), Ad_g1inv
   disp('Found  '), Ad_g1_inv
else
   disp('Test 2 OK.')
end

%% A second twist transformed by the first
w2 = [0;0;1];
v2 = [2;0;0];
tw2 = [v2;w2];

g1 = expr(tw1, pi/2);

[twn2, Ad_g1] = adjoint_trf(tw2, g1);

if norm(twn2(4:6) - tw2(4:6)) > thr
   disp('Test 3 failed.')
   disp('Expected  '), tw2
   disp('Found  '), twn2
else
   disp('Test 3 OK.')
end

if norm(twn2(1:3) - [1;1;0]) > thr
   disp('Test 4 failed.')
   disp('Expected  '), [1;1;0]
   disp('Found  '), twn2(1:3)
   keyboard
else
   disp('Test 4 OK.')
end

