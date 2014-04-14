function visualize_mobility_FT(WW, colors, clubpath, fig, fname)
%%  visualize_mobility(W, colors, fig, fname)
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

figure(fig)
clf

set(fig, 'position', [100, 100, 700, 350])
nmob = size(WW,3);
for i = 1:nmob
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

    axlim1 = [-4, 4];
    axlim2 = [-2, 8];
    
    subplot(1,3,1)
    plot(eexz(1,:), eexz(2,:), 'color', color, 'linewidth', 3)
    hold on
    plot(3*clubpath(:,1), 3*clubpath(:,3), ':ob')
    title('Sagittal plane')
    xlabel('x')
    ylabel('z')
    axis equal
    set(gca, 'xlim', axlim1)
    set(gca, 'ylim', axlim2)
    annotation('ellipse', [0.2262 0.2537 0.017 0.02], 'EdgeColor',[0 0 0]);
    annotation('ellipse', [0.2262+0.006 0.2537+0.008 0.0049 0.005], 'FaceColor',[0 0 0]);
    axis off
        
    subplot(1,3,2)
    plot(eeyz(1,:), eeyz(2,:), 'color', color, 'linewidth', 3)
    hold on
    plot(3*clubpath(:,2), 3*clubpath(:,3), ':ob')
    title('Frontal plane')
    xlabel('y')
    ylabel('z')
    axis equal
    set(gca, 'xlim', axlim1)
    set(gca, 'ylim', axlim2)
   annotation('arrow',[0.4578 0.5776], [0.2637 0.2637]);
    axis off
        
    subplot(1,3,3)
    plot(eexy(2,:), -eexy(1,:), 'color', color, 'linewidth', 3)
    hold on
    plot(3*clubpath(:,2), -3*clubpath(:,1), ':ob')
    title('Horizontal plane')
    xlabel('y')
    ylabel('x')
    axis equal
    set(gca, 'xlim', axlim1)
    set(gca, 'ylim', axlim2)
    annotation('arrow',[0.7384 0.8582], [0.2637 0.2637]);
    axis off    
end

print(fname, '-dpdf')

disp(['Wrote figure to file ', fname])
