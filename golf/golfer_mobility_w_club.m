function [W, Wep, M, Mep] = golfer_mobility_w_club(bm, bstates, cm, cstates)
%%   [W, Wep, M, Mep] = golfer_mobility_w_club(bm, bstates, cm, cstates)
%% Returns the mobility of the golfer end segment (whole club), a (6 x 6 x nfrs) matrix
%%  and the mobility of the end point (club head).
%%
%% Input
%%    bm          ->  kinematic model of body
%%    bstates      ->  (nst x nfrs) sequence of joint angles.
%%    cm          ->  kinematic model of club
%%    cstates      ->  (6 x nfrs) generalized coordinates for club
%% Output
%%    W           <-  the whole club mobility (6 x 6 x nfrs),       
%%    Wep         <-  the mobility of the club head (3 x 3 x nfrs),       

%% Kjartan Halvorsen
%% 2013-09-13

if (nargin == 0)
  do_unit_test();
  return
end

[nst, nfrs] = size(bstates);

objdof = 6;
epdof = 3;

%% Get the link jacobian
Jh = end_link_jacobian(bm, bstates); % ( (6+6) x nst x nfrs ) jacobian of end effector body vel
Jh_ep = end_point_jacobian(bm, bstates); % ( (3 + 3) x nst x nfrs ) jacobian in spatial vel 
Jhc = end_link_jacobian(cm, cstates); % ( 6 x nst x nfrs ) jacobian of end effector body vel
Jhc_ep = end_point_jacobian(cm, cstates); % ( 3 x nst x nfrs ) jacobian in spatial vel 
 
%% Get the link inertia
Mf = generalized_manipulator_inertia(bm, bstates); % (nst x nst x nfrs)
Mfc = generalized_manipulator_inertia(cm, cstates); % (nst x nst x nfrs)

%% Form Jbar and Gbar
%% Check for two or one endlink.
nendlinks = size(Jh,1) / 6;

if nendlinks == 2
  Gtransp = repmat(eye(6), [2, 1]);
  Gtransp_ep = repmat(eye(3), [2, 1]);
else
  Gtransp = eye(6);
  Gtransp_ep = eye(3);
end

%keyboard
%% Compute the mass matrix and invert it
W = zeros(6, 6, nfrs);
Wep = zeros(3, 3, nfrs);
M = zeros(6, 6, nfrs);
Mep = zeros(3, 3, nfrs);

for i = 1:nfrs
    nsteff = rank(Jh(:,:,i)');
    dofdif = 12-nsteff;
    if nendlinks == 2
      Gbartransp = [Gtransp  zeros(12, nst-2*objdof+dofdif)
		    zeros(nst-2*objdof+dofdif,6) eye(nst-2*objdof+dofdif, nst-2*objdof+dofdif)];
      nsteff = rank(Jh_ep(:,:,i)');
      dofdif = 6-nsteff;
      Gbartransp_ep = [Gtransp_ep  zeros(6, nst-2*epdof+dofdif)
		       zeros(nst-2*epdof+dofdif,3) eye(nst-2*epdof+dofdif, nst-2*epdof+dofdif)];
    else
	Gbartransp = [Gtransp zeros(6, nst-objdof)
		      zeros(nst-objdof, 6) eye(nst-objdof)];
	Gbartransp_ep = [Gtransp_ep zeros(3, nst-3)
			 zeros(nst-3, 3) eye(nst-3)];
    end
    Jbar = cat(1, Jh(:,:,i), (null(Jh(:,:,i)))' );
    Jbar_ep = cat(1, Jh_ep(:,:,i), (null(Jh_ep(:,:,i)))' );
    Jbarinv = pinv(Jbar);
    Jbarcinv = pinv(Jhc(:,:,i));
    %%keyboard
    Mbar = Gbartransp'*Jbarinv'*Mf(:,:,i)*Jbarinv*Gbartransp;
    Mbar(1:6, 1:6) = Mbar(1:6, 1:6)  + Jbarcinv'*Mfc(:,:,i)*Jbarcinv;
    M(:,:,i) = Mbar(1:6, 1:6);
    WW = pinv(Mbar);
    %WW = pinv(Mbar(1:6, 1:6));
    W(:,:,i) = WW(1:6, 1:6);

    Jbarinv = pinv(Jbar_ep);
    Jbarcinv = pinv(Jhc_ep(:,:,i));
    Mbar = Gbartransp_ep'*Jbarinv'*Mf(:,:,i)*Jbarinv*Gbartransp_ep;
    Mbar(1:3, 1:3) = Mbar(1:3, 1:3) + Jbarcinv'*Mfc(:,:,i)*Jbarcinv;
    Mep(:,:,i) = Mbar(1:3,1:3);
    WW = pinv(Mbar);
    %%WW = pinv(Mbar(1:3, 1:3));
    Wep(:,:,i) = WW(1:3, 1:3);
end

function do_unit_test()

  %% Load some data and generate model
  %% Two planar arms
  pos1 = cat(2, ...
	     [-2;0;0],...
	     [-2;1;0],...
	     [-2;2;0],...
	     [-1;2;0],...
	     [0;2;0]);
  pos2 = cat(2, ...
	     [2;0;0],...
	     [2;1;0],...
	     [2;2;0],...
	     [1;2;0],...
	     [0;2;0]);

  pos1 = cat(2, ...
	     [-1;0;0],...
	     [-1;1;0],...
	     [-1;2;0],...
	     [-1;3;0]);

  pos2 = cat(2, ...
	     [1;0;0],...
	     [1;1;0],...
	     [1;2;0],...
	     [1;3;0]);

  masses = ones(4,1);

  g0 = eye(4);
  object_frame = [eye(3) [0;3;0]; zeros(1,3) 1];
  
  
  pm1 = link_model(pos1, masses, [], g0, object_frame);
  pm2 = link_model(pos2, masses, [], g0, object_frame);

  pm = combine_models([], pm1, pm2);

  states = zeros(18,1);
  %states([4 8]) = 0.3;

  isplanar = 1;
  [W, Wep, M, Mep] = golfer_mobility(pm, states);
  keyboard
  
  
  
  
