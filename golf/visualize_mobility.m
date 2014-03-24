function visualize_mobility(WW, colors, fig, fname)
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

nmob = size(WW,3);
for i = 1:nmob
    W = WW(:,:,i);
    color = colors(i,:);

    Wxz = W([1 3], [1 3]);
    Wyz = W([2 3], [2 3]);
    [Vxz,Dxz] = eig(Wxz);
    [Vyz,Dyz] = eig(Wyz);
    VVxz = Vxz*sqrt(Dxz);
    VVyz = Vyz*sqrt(Dyz);

    t = linspace(0, 2*pi, 400);
    ee = [cos(t); sin(t)];  % Unit circle

    
    eexz = VVxz*ee;
    eeyz = VVyz*ee;


    axlim = [-2.5, 2.5];

    subplot(121)
    plot(eexz(1,:), eexz(2,:), 'color', color, 'linewidth', 3)
    hold on
    title('xz-plane')
    set(gca, 'xlim', axlim)
    set(gca, 'ylim', axlim)
    xlabel('x')
    ylabel('z')
    axis equal
    
    subplot(122)
    plot(eeyz(1,:), eeyz(2,:), 'color', color, 'linewidth', 3)
    hold on
    title('yz-plane')
    set(gca, 'xlim', axlim)
    set(gca, 'ylim', axlim)
    xlabel('y')
    ylabel('z')
    axis equal
end

print(fname, '-dpdf')

disp(['Wrote figure to file ', fname])
