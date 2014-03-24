function [twn, Ad_g_inv] = adjoint_trf_inv(tw, g)
%%  [twn, Ad_g_inv] = adjoint_trf_inv(tw, g)
%% Transforms the twist tw into twn by applying the inverse adjoint transformation associated with g.
%% See eq 2.58 in Murray, Li, Sastry:
%%   Ad_g = [R'  -\hat{R'p}R' ]
%%          [0      R'    ],
%% where 
%%   g = [R  p]
%%       [0  1].

%% Kjartan Halvorsen
%% 2013-05-31

  if size(tw,1) == 4
    twv = vee(tw);
  else
    twv = tw;
  end

  R = g(1:3, 1:3);
  p = g(1:3, 4);

  Ad_g_inv = [ R'  -hat(R'*p)*R'
	       zeros(3,3) R'];

  twn = Ad_g_inv*twv;

  if size(tw, 1) == 4
    twn = hat(twn);
  end
end

