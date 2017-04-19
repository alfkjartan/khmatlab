function visualize_mobility_FT_doubleaxes(WW, colors, clubpath, figs, fnames)
%%  Based on visualize_mobility_FT(W, colors, figs, fnames), but
%% NO SCALING OF WW OR PATH!! This is handled by the double axes in the same plots
%%
%%  Plots the mobility as an ellipse in the xz-plane (frontal plane of golfer) and the yz-plane
%%  (sagittal plane of golfer)
%%
%%  Input
%%     W      ->  The mobility matrix at impact. If third dimension, nmob > 1, plot several ellipses
%%     colors ->  Specify colors for the mobility matrices (nmob x 3).
%%     fig    ->  The figure to plot in
%%     fname  ->  The name to write pdf to

%% Kjartan Halvorsen
%% 2014-02-04

fszax = 14;
fsztxt = 16;
for i=1:length(figs)
  figure(figs(i))
  clf
end

nmob = size(WW,3);
haxes2s = [];
haxes1s = [];
for i = 1:nmob
    
    figure(figs(1))
    set(figs(1), 'position', [100, 100, 350, 350])
    W = WW(:,:,i);
    color = colors(i,:);

    Wxz = W([1 3], [1 3]);
    Wxy = W([1 2], [1 2]);
    Wyz = W([2 3], [2 3]);
    [Vxz,Dxz] = eig(Wxz);
    [Vxy,Dxy] = eig(Wxy);
    [Vyz,Dyz] = eig(Wyz);
    VVxz = Vxz*sqrt(Dxz);
    VVxy = Vxy*sqrt(Dxy);
    VVyz = Vyz*sqrt(Dyz);

    t = linspace(0, 2*pi, 400);
    ee = [cos(t); sin(t)];  % Unit circle

    
    eexz = VVxz*ee;
    eexy = VVxy*ee;
    eeyz = VVyz*ee;
    %keyboard

    axlim2_x = [-4, 4];
    axlim2_y = [-4, 4];
    axlim1_x = [-1, 1];
    axlim1_y = [-1, 1];
    
    if length(haxes1s) < 1
      line(clubpath(:,1), clubpath(:,3), 'LineStyle', 'none', 'Marker', 'o')
      axis equal
      haxes1 = gca;
      haxes1s = cat(1, haxes1s, haxes1);
      haxes1_pos = get(haxes1,'Position') % store position of first axes
      haxes1_pos(2) = 0.75*haxes1_pos(2);
      set(haxes1, 'Position', haxes1_pos);
      set(haxes1, 'FontSize', fszax)
      set(haxes1, 'xlim', axlim1_x)
      set(haxes1, 'ylim', axlim1_y)
      xlabel('x [m]')
      ylabel('z [m]')
      hold on
    else
	haxes1 = haxes1s(1);
    end

    if length(haxes2s) < 1
      haxes1_pos = get(haxes1,'Position') % store position of first axes
      haxes2 = axes('Position',haxes1_pos,...
		    'XAxisLocation','top',...
		    'YAxisLocation','right',...
		    'Color','none', 'FontSize', fszax);
      haxes2s = cat(1, haxes2s, haxes2);
    else
	haxes2 = haxes2s(1);
    end
    line(eexz(1,:), eexz(2,:), 'Parent', haxes2,'color', color, 'linewidth', 3)
    %title('Sagittal plane')
    hold on
	
    xlabel('x [kg^{-1}]')
    ylabel('z [kg^{-1}]')
    axis equal
    set(haxes2, 'xlim', axlim2_x)
    set(haxes2, 'ylim', axlim2_y)

    %annotation('ellipse', [0.2272 0.2237 0.02 0.028], 'EdgeColor',[0 0 0]);
    %h=text(-0.2, -4.9, 'X','Color',[0 0 0]);
    %get(h)
    %axis off
        
    figure(figs(2))
    set(figs(2), 'position', [200, 550, 350, 350])

    if length(haxes1s) < 2
      line(clubpath(:,2), clubpath(:,3), 'LineStyle', 'none', 'Marker', 'o')
      axis equal
      haxes1 = gca;
      haxes1s = cat(1, haxes1s, haxes1);
      haxes1_pos = get(haxes1,'Position') % store position of first axes
      haxes1_pos(2) = 0.75*haxes1_pos(2);
      set(haxes1, 'Position', haxes1_pos);
      set(haxes1, 'FontSize', fszax)
      set(haxes1, 'xlim', axlim1_x)
      set(haxes1, 'ylim', axlim1_y)
      xlabel('y [m]')
      ylabel('z [m]')
      hold on
    else
      haxes1 = haxes1s(2);
    end


    if length(haxes2s) < 2
      haxes1_pos = get(haxes1,'Position'); % store position of first axes
      haxes2 = axes('Position',haxes1_pos,...
		    'XAxisLocation','top',...
		    'YAxisLocation','right',...
		    'Color','none', 'FontSize', fszax);
      haxes2s = cat(1, haxes2s, haxes2);
    else
	haxes2 = haxes2s(2);
    end
    line(eeyz(1,:), eeyz(2,:), 'Parent', haxes2,'color', color, 'linewidth', 3)
    hold on
    %title('Frontal plane')
    xlabel('y [kg^{-1}]')
    ylabel('z [kg^{-1}]')
    axis equal
    set(haxes2, 'xlim', axlim2_x)
    set(haxes2, 'ylim', axlim2_y)

    %annotation('arrow',[0.4578 0.5776], [0.2237 0.2237]);
    %axis off
        
    figure(figs(3))
    set(figs(3), 'position', [250, 900, 350, 350])

    if length(haxes1s) < 3
      line(clubpath(:,2), clubpath(:,1), 'LineStyle', 'none', 'Marker', 'o')
      axis equal
      haxes1 = gca;
      haxes1s = cat(1, haxes1s, haxes1);
      haxes1_pos = get(haxes1,'Position') % store position of first axes
      haxes1_pos(2) = 0.75*haxes1_pos(2);
      set(haxes1, 'Position', haxes1_pos);
      set(haxes1, 'FontSize', fszax)
      set(haxes1, 'xlim', axlim1_x)
      set(haxes1, 'ylim', axlim1_y)
      set(haxes1,'YDir','reverse');
      xlabel('y [m]')
      ylabel('x [m]')
      hold on
    else
      haxes1 = haxes1s(3);
    end

    if length(haxes2s) < 3
      haxes1_pos = get(haxes1,'Position'); % store position of first axes
      haxes2 = axes('Position',haxes1_pos,...
		    'XAxisLocation','top',...
		    'YAxisLocation','right',...
		    'Color','none', 'FontSize', fszax);
      haxes2s = cat(1, haxes2s, haxes2);
    else
	haxes2 = haxes2s(3);
    end
    line(eexy(1,:), eexy(2,:), 'Parent', haxes2,'color', color, 'linewidth', 3)
    set(haxes2,'YDir','reverse');
    hold on
    %title('Horizontal plane')
    xlabel('y [kg^{-1}]')
    ylabel('x [kg^{-1}]')
    axis equal
    set(haxes2, 'xlim', axlim2_x)
    set(haxes2, 'ylim', axlim2_y)

    %annotation('arrow',[0.7384 0.8582], [0.2237 0.2237]);
    %axis off    
end

for i=1:length(figs)
  figure(figs(i))
  set(findall(gcf,'type','text'),'FontSize',fsztxt)
  print(fnames{i}, '-dpdf')

  disp(['Wrote figure to file ', fnames(i)])
end
