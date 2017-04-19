     if debug
       % Plot mobility (Wep) in the three directions
	figure(mobfig(1));
	clf
	%Wepx = Wep(1,1,:);

	lstfr = min(400, size(Wepbothleft,3));
	subplot(311)
	plot(vectorize(Wepbothleft(1,1,1:lstfr)), 'color', [1, 0, 0]);
	hold on
	plot(vectorize(Wepbothright(1,1,1:lstfr)), 'color', [0 1 0]);
	subplot(312)
	plot(vectorize(Wepbothleft(2,2,1:lstfr)), 'color', [1, 0, 0]);
	hold on
	plot(vectorize(Wepbothright(2,2,1:lstfr)), 'color', [0 1 0]);
	subplot(313)
	plot(vectorize(Wepbothleft(3,3,1:lstfr)), 'color', [ 1, 0, 0 ]);
	hold on
	plot(vectorize(Wepbothright(3,3,1:lstfr)), 'color', [0 1 0]);

	figure(mobfig(2));
	clf
	%Wepx = Wep(1,1,:);

	lstfr = min(410, size(Wbothleft,3));
	subplot(311)
	plot(vectorize(Wbothleft(1,1,1:lstfr)), 'color', [1, 0, 0]);
	hold on
	plot(vectorize(Wbothright(1,1,1:lstfr)), 'color', [0 1 0]);
	subplot(312)
	plot(vectorize(Wbothleft(2,2,1:lstfr)), 'color', [1, 0, 0]);
	hold on
	plot(vectorize(Wbothright(2,2,1:lstfr)), 'color', [0 1 0]);
	subplot(313)
	plot(vectorize(Wbothleft(3,3,1:lstfr)), 'color', [ 1, 0, 0 ]);
	hold on
	plot(vectorize(Wbothright(3,3,1:lstfr)), 'color', [0 1 0]);

	figure(inertfig);
	clf
	%Wepx = Wep(1,1,:);

	subplot(311)
	plot(vectorize(Wepleft(1,1,1:lstfr)), 'color', [1, 0, 0]);
	hold on
	plot(vectorize(Wepright(1,1,1:lstfr)), 'color', [0 1 0]);
	plot(vectorize(Wepbothright(1,1,1:lstfr)), 'color', [0 0 0]);
	subplot(312)
	plot(vectorize(Wepleft(2,2,1:lstfr)), 'color', [1, 0, 0]);
	hold on
	plot(vectorize(Wepright(2,2,1:lstfr)), 'color', [0 1 0]);
	plot(vectorize(Wepbothright(2,2,1:lstfr)), 'color', [0 0 0]);
	subplot(313)
	plot(vectorize(Wepleft(3,3,1:lstfr)), 'color', [ 1, 0, 0 ]);
	hold on
	plot(vectorize(Wepright(3,3,1:lstfr)), 'color', [0 1 0]);
	plot(vectorize(Wepbothright(3,3,1:lstfr)), 'color', [0 0 0]);

     end

    % Find events. 
     [imp_fr, imp_fit, back_starts, back_ends] = find_events_new(objdboth(:,1:3));
     
    if ~isempty(imp_fit)
            imp_fr = round(imp_fit);
    end
    
    % Find events. 
    %events{1,1} = 'impact';
    % events{1,2} = imp_fr;
    % events{2,1} = 'backswingstarts';
    % events{2,2} = back_starts;
    % events{3,1} = 'backswingends';
    % events{3,2} = back_ends;

    % Does nothing right now. Not implemented
    Wepboth = Wepbothleft;
    [Wpathleft, Wnormleft] = mobility_along_path(Wepleft, objdleft);
    [Wpathright, Wnormright] = mobility_along_path(Wepright, objdright);
    [Wpathboth, Wnormboth] = mobility_along_path(Wepboth, objdboth(:,1:3));

    mobility_along_path_at_impact_left(fpind, trind) = Wpathleft(imp_fr);
    mobility_normal_to_path_at_impact_left = Wnormleft(imp_fr);

    mobility_along_path_at_impact_right(fpind, trind) = Wpathright(imp_fr);
    mobility_normal_to_path_at_impact_right = Wnormright(imp_fr);

    mobility_along_path_at_impact_both(fpind, trind) = Wpathboth(imp_fr);
    mobility_normal_to_path_at_impact_both = Wnormboth(imp_fr);


    [pth,mfname] = fileparts(filestr);

    fname_plot = fullfile(ndir, [mfname,'_mobility.pdf']);
    fname_plot_imp = fullfile(ndir, [mfname,'_mobility_impact.pdf']);


     if plot_mobility
	%% Ellipse at impact
	Watimpact = cat(3, Wepboth(:,:,imp_fr), Wepleft(:,:,imp_fr), Wepright(:,:,imp_fr));
	colors = [0 0 0; 0.9 0 0; 0 0.8 0];
	visualize_mobility(Watimpact, colors, mobfig(3), fname_plot_imp);

	
       % Plot mobility (Wep) in the three directions
	strtfr = back_starts;
	endfr = imp_fr + 20;

	impact = imp_fr - strtfr + 1;
	figure(mobfig(2));
	clf

	subplot(211)
	ylabel('W along path')
	plot(Wpathleft(strtfr:endfr), 'color', [0.9, 0, 0]);
	hold on
	plot(Wpathright(strtfr:endfr), 'color', [0, 0.8, 0]);
	plot(Wpathboth(strtfr:endfr), 'color', [0, 0, 0]);
	yl = get(gca, 'Ylim');
	plot([impact, impact], yl, 'color', [0 0 0])

	text(impact+5, 0.95*yl(2), sprintf('Left: %1.2f', Wpathleft(imp_fr)), 'color', [0.9, 0, 0])
	text(impact+5, 0.9*yl(2), sprintf('Right: %1.2f', Wpathright(imp_fr)), 'color', [0, 0.8, 0])
	text(impact+5, 0.85*yl(2), sprintf('Both: %1.2f', Wpathboth(imp_fr)), 'color', [0.7, 0.7, 0])

	subplot(212)
	ylabel('W normal to path')
	plot(Wnormleft(strtfr:endfr), 'color', [0.9, 0, 0]);
	hold on
	plot(Wnormright(strtfr:endfr), 'color', [0, 0.8, 0]);
	plot(Wnormboth(strtfr:endfr), 'color', [0, 0, 0]);
	yl = get(gca, 'Ylim');
	plot(repmat(imp_fr-strtfr+1, 1, 2), yl, 'color', [0 0 0])

	text(impact+5, 0.95*yl(2), sprintf('Left: %1.2f', Wnormleft(imp_fr)), 'color', [0.9, 0, 0])
	text(impact+5, 0.9*yl(2), sprintf('Right: %1.2f', Wnormright(imp_fr)), 'color', [0, 0.8, 0])
	text(impact+5, 0.85*yl(2), sprintf('Both: %1.2f', Wnormboth(imp_fr)), 'color', [0, 0, 0])

	print(fname_plot, '-dpdf')



     end
   
