function [imp_fr, imp_fit, backswingstarts, backswingends] = find_events_new(clubheadcenter)
%  [imp_fr, imp_fit, backswingstart, backswingends]  = find_events_new(objectcetner)
% Function that finds the impact event given the golf model gm, and
% the estimated states.
% 
% Input
%    objectcenter  ->  Trajectort of club head center (nfrs x 3)
% Output
%    imp_fr        <-  The impact event, in frame number
%    imp_fit       <-  The impact event, interpolated
%    backswingstarts  <-  The start of the backswing
%    backswingends    <-  The end of the backswing (top)

% Kjartan Halvorsen
% 2009-07-11   Eirik fyller 3år!!!!

% Revisions
% 2014-02-03    Assume path of clubhead center already computed elsewhere, 
%               and passed as argument.

useOld = 1;
askForImpact = 0;

debug = 0;

% The interesting directions
e_goal = [0;1;0];
e_ap   = [1;0;0];
e_vert = [0;0;1];

% Get the movement of the clubhead

ep = clubheadcenter;
% The endpoint path and mechanism Jacobian
nfrs = size(ep,1);

ep_vel = centraldiff(ep, 120); % The endpoint velocity

if useOld
  if askForImpact
    inputanswers = ...
	inputdlg({['Enter x-coordinate for club toe at impact: '], ...
		  'Enter the first frame number: '});
    
    
    impactycoord = str2num(inputanswers{1});
    firstframe = str2num(inputanswers{2});
  else
    impactycoord = ep(1,2);
  end

  [impact, impact_fit, pquad, dist2address, firstmax, backswingstarts] = find_impact_from_point_new(ep);
  imp_fr = impact;
  imp_fit = impact_fit;


  clubtoe = ep;
  
  clubtoevelocity  =  centraldiff(ep, 120); % The endpoint velocity

  %keyboard
  
  magnvel = sqrt(sum(clubtoevelocity.^2, 2));
  
  % End of backswing
  clubhighbeforeimpact = intersect(find(clubtoe(:,3) > -0.9), (1:impact)');
  clubvelmin = min(magnvel(clubhighbeforeimpact));
  backswingends = find(magnvel == clubvelmin(1));
    
else

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

end

if debug
  
  figure(40)
  clf
  xstart = imp_fr-4;
  xend = imp_fr+4;
  xplot = xstart:xend;
  xplot2 = linspace(xstart,xend,100);
  fittedcurve = polyval(pcub, xplot2-imp_fr);
  plot(xplot, dist2address(xplot));
  hold on
  plot(xplot2, fittedcurve, 'r');
  %keyboard
  plot([imp_fit imp_fit], get(gca,'YLim'), 'g')
end

  



