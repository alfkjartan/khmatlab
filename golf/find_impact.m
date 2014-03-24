function [imp_fr, imp_fit] = find_impact(gm, states)
%  imp = find_impact(gm, states)
% Function that finds the impact event given the golf model gm, and
% the estimated states.
% 
% Input
%    gm       ->  golf model, as returned by build_golf_model.m
%    states   ->  (nst x nfr) sequence of states.
% Output
%    imp      <-  The impact event

% Kjartan Halvorsen
% 2009-07-11   Eirik fyller 3år!!!!

debug = 1;

% The interesting directions
e_goal = [0;1;0];
e_ap   = [1;0;0];
e_vert = [0;0;1];

% Get the movement of the clubhead

% The endpoint path and mechanism Jacobian
nst = size(states,1)/2;
nfrs = size(states,2);

[slask1, slask2, ep, epnames] = sim_model(gm,states(1:nst,:), 'CoM');
% The above function will return two end points, due to the fact
% that the club is at the end of two chains. Use only one
ep = ep(:,1:3);

ep_vel = centraldiff(ep, 120); % The endpoint velocity

% Find the club position at address. To do this, first find the
% take away event.
threshold = 0.2; % m/s
frames_after = 20;
tw = threshold_event(ep_vel, -e_goal, threshold, frames_after);

% Club head position nearest ball is the maximal position in the
% direction of the goal, prior to take away.
[address_pos, address_pos_frame] = find_max(ep, e_goal, 1, tw);

% Compute distance from current club head position to address_pos
dist2address = repmat(address_pos, nfrs, 1) - ep;
dist2address = sqrt(sum(dist2address.^2,2));

% Find the event when this distance is minimal.
%imp_temp = tw + ...
%    find(dist2address(tw+1:end) == min(dist2address(tw+1:end)));
[dist, imp_temp] = find_max(-dist2address, [], tw);

% Interpolate by cubic spline using  two frames before and two after
xt = (imp_temp-2):(imp_temp+2);
xtc = -2:2;
pcub = polyfit(xtc', dist2address(xt), 3);
pderiv = [3*pcub(1) 2*pcub(2) pcub(3)];
proots = roots(pderiv);

%large_roots = find(proots > imp_temp-1);
%small_roots = find(proots < imp_temp+1);
large_roots = find(proots > -1);
small_roots = find(proots < +1);

the_root_ind = intersect(large_roots, small_roots);
imp_fit = imp_temp + proots(the_root_ind);
imp_fr = imp_temp;

if debug
  
  figure(40)
  clf
  xstart = imp_temp-4;
  xend = imp_temp+4;
  xplot = xstart:xend;
  xplot2 = linspace(xstart,xend,100);
  fittedcurve = polyval(pcub, xplot2-imp_temp);
  plot(xplot, dist2address(xplot));
  hold on
  plot(xplot2, fittedcurve, 'r');
  plot([imp imp], get(gca,'YLim'), 'g')
end

  



