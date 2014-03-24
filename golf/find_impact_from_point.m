function  [impact, impact_fit, pquad, dist2address, firstmax, backswingstarts] = find_impact_from_point(ep)

				% Find events
  clubtoe = ep;
  
  clubtoevelocity  =  centraldiff(ep, 120); % The endpoint velocity

  %keyboard

  magnvel = sqrt(sum(clubtoevelocity.^2, 2));
  
  clubtoehasdata = find(clubtoe(:,2) ~= 0);
  
				% Start of backswing
  clublow = intersect(clubtoehasdata, find(clubtoe(:,3) < -0.9));

  clublowandslow = intersect(... 
			     intersect(clublow, find(clubtoevelocity(:,2) > -1)), ...
			     find(magnvel < 1));
  backswingstarts = clublowandslow(end)-1;
				% Impact
  % redefine impactycoord as max ycoord prior to start of backswing
  impactycoord = max(clubtoe(1:backswingstarts, 2));

  % The first maximum before backswingstarts is also of interest
  firstmax = intersect((1:backswingstarts), ...
		       find( clubtoevelocity(:,2) < 0 ));
  firstmax = firstmax(end);


  clubback = find(clubtoe(:,2) < (impactycoord - 0.1));

  clubbackandlow = intersect(clublow, clubback);
  clubhighvel = find(magnvel > 10);
  club_low_and_highvel = intersect(clublow, clubhighvel);
  clubbefore = intersect(club_low_and_highvel, ...
		       find(clubtoe(:,2) < impactycoord));
  clubafter = intersect(club_low_and_highvel, ...
		      find(clubtoe(:,2) > ...
			   impactycoord));
%keyboard
  impact_frs = [clubbefore(end) clubafter(1)]
  if (abs(clubtoe(impact_frs(1), 1) - impactycoord) < ...
      abs(clubtoe(impact_frs(2), 1) - impactycoord))
    impact = impact_frs(1);
  else
    impact = impact_frs(2);
  end

  % Interpolate by cubic spline using  four frames before
  xt = ((impact_frs(1)-3):impact_frs(1));
  xtc = -3:0;
  dist2address = impactycoord - clubtoe(:,2);
  pquad = polyfit(xtc', dist2address(xt), 3);
  %pderiv = [3*pcub(1) 2*pcub(2) pcub(3)];
  proots = roots(pquad);

  %large_roots = find(proots > imp_temp-1);
  %small_roots = find(proots < imp_temp+1);
  large_roots = find(proots > -1);
  small_roots = find(proots < +1);

  the_root_ind = intersect(large_roots, small_roots);
  impact_fit = impact_frs(1) + proots(the_root_ind);
