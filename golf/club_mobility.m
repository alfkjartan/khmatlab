function [W, Wep] = club_mobility(gm, states)
%%  [W, Wep] = club_mobility(gm, states)
%% Returns the mobility of the whole club (6 x 6 x nfrs) matrix and the mobility of the 
%% end point (club head).
%%
%% Input
%%    gm          ->  golf model, constructed by build_golf_model_w_inertia
%%    states      ->  (nst x nfrs) sequence of joint angles.
%% Output
%%    W           <-  the whole club mobility (6 x 6 x nfrs),       
%%    Wep         <-  the mobility of the club head (3 x 3 x nfrs),       

%% Kjartan Halvorsen
%% 2013-06-13

%% Get the link jacobian
Jh = link_jacobian(gm, states); % ((6+6) x nst x nfrs)

%% Get the link inertia
Mf = generelized_manipulator_inertia(gm, states); % (nst x nst x nfrs)

%% For Jbar
Jbar = 
