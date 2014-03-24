
%% PANEL DEMO 2
%%
%% (a) Create an empty panel object
%% (b) Packing multiple panels, using different edges
%% (c) Plotting into those panels, setting axis labels
%% (d) Control render metrics, for a specific panel



%% (a)

% clear or create figure and create root panel covering
% whole figure
clf
p = panel;



%% (b)

% form an upper panel A and a lower panel B where A takes 30% of the
% available space, leaving 70% for B
pA = p.pack(30);
pB = p.pack();

% form sub-panels within B
pB.edge = 'l'; % pack on the left edge (default is to pack on the top edge)
pBA = pB.pack(25); % use leftmost 25% of the space for BA
pBB = pB.pack(25); % use next 25% of the space for BB
pBC = pB.pack(); % use remaining for BC
pBBA = pBC.pack(1/3); % values less than 1 are assumed to be fractions rather than percentages
pBBB = pBC.pack(1/3); % you don't have to use the 'remainder' space
pBBC = pBC.pack(1/3); % here, we've used all the space explicitly anyway



%% (c)

% plot into panel
select(pA)
x = -100:100;
plot(x,x.^2/1000,'r')
set(gca,'xtick',[-75:25:75])
pA.xlabel = 'x is great';
pA.ylabel = 'y is better';

% plot into panel
select(pBA)
plot(randn(1000,3))
pBA.xlabel = 'foo';
pBA.ylabel = 'bar';
x = 1:10;

% plot into panel
select(pBB)
plot(randn(1000,3))
pBB.xlabel = 'foo';
pBB.ylabel = 'bar';
x = 1:10;

% plot into panel
select(pBBA)
plot(x,x)
pBBA.xlabel = 'x';
pBBA.ylabel = 'x';

% plot into panel
select(pBBB)
plot(x,x.^2)
pBBB.xlabel = 'x';
pBBB.ylabel = 'x^2';

% plot into panel
select(pBBC)
plot(x,x.^3)
pBBC.xlabel = 'x';
pBBC.ylabel = 'x^3';



%% (d)

% need some extra margin (on my monitor, anyway!)
% - try it without, see what's different!
pBC.parentmargin = [8 0 0 0]; % panel interface is in mm, by default

