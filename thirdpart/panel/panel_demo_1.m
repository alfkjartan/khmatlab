
%% PANEL DEMO 1
%%
%% (a) Create an empty panel object
%% (b) Packing two panels into one figure
%% (c) Plotting into those panels, setting axis labels
%% (d) Control render metrics, figure-wide



%% (a)

% clear or create figure and create root panel covering
% whole figure
clf
p = panel;



%% (b)

% form an upper panel A and a lower panel B, which will
% share the unassigned space in p.
pA = p.pack();
pB = p.pack();



%% (c)

% plot into panel A
select(pA)
x = -100:100;
plot(x,x.^2/1000,'r')

% set axis labels through panel, for control over metrics
pA.xlabel = 'x is great';
pA.ylabel = 'y is better';

% plot into panel B
select(pB)
plot(randn(1000,3))
axis([0 1000 -5 5])

% set axis labels through panel, for control over metrics
pB.xlabel = 'foo';
pB.ylabel = 'bar';



%% (d)

% set render properties on root panel - properties are
% inherited by all children and grandchildren, etc.
p.fontname = 'times';
p.fontsize = 10;
p.fontweight = 'bold';

