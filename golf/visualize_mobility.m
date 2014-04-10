function visualize_mobility(WW, colors, , endpointpath, fig, fname)
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

nmob = size(WW,3);
endpointpathscale = 6;
if nargin < 5
   endpointpath = repmat(nan, [2, 3, nmob]);
end

figure(fig)
clf

lwdth = 5;

for i = 1:nmob
    W = WW(:,:,i);
    color = colors(i,:);

    Wxz = W([1 3], [1 3]);
    Wyz = W([2 3], [2 3]);
    Wxy = W([1 2], [1 2]);
    [Vxz,Dxz] = eig(Wxz);
    [Vyz,Dyz] = eig(Wyz);
    [Vxy,Dxy] = eig(Wxy);
    VVxz = Vxz*sqrt(Dxz)/endpointpathscale;
    VVyz = Vyz*sqrt(Dyz)/endpointpathscale;
    VVxy = Vxy*sqrt(Dxy)/endpointpathscale;

    t = linspace(0, 2*pi, 400);
    ee = [cos(t); sin(t)];  % Unit circle

    
    eexz = VVxz*ee;
    eeyz = VVyz*ee;
    eexy = VVxy*ee;


    axlim = [-1, 1];

    subplot(221)
    plot(eeyz(1,:), eeyz(2,:), 'color', color, 'linewidth', lwdth)
    hold on
    plot(endpointpath(:,2,i), endpointpath(:,3,i), 'color', color, 'linewidth', lwdth-1)
    title('Frontal plane (yz)')
    set(gca, 'xlim', axlim)
    set(gca, 'ylim', axlim)
    xlabel('y')
    ylabel('z')

    axis equal
    
    subplot(222)
    plot(eexz(1,:), eexz(2,:), 'color', color, 'linewidth', lwdth)
    hold on
    plot(endpointpath(:,1,i), endpointpath(:,3,i), 'color', color, 'linewidth', lwdth-1)
    title('Sagittal plane (xz)')
    set(gca, 'xlim', axlim)
    set(gca, 'ylim', axlim)
    xlabel('x')
    ylabel('z')
    axis equal

    subplot(223)
    plot(eexy(1,:), eexy(2,:), 'color', color, 'linewidth', lwdth)
    hold on
    plot(endpointpath(:,1,i), endpointpath(:,2,i), 'color', color, 'linewidth', lwdth-1)
    title('Horizontal plane (xy)')
    set(gca, 'xlim', axlim)
    set(gca, 'ylim', axlim)
    xlabel('x')
    ylabel('y')
    axis equal
end

print(fname, '-dpdf')

disp(['Wrote figure to file ', fname])
