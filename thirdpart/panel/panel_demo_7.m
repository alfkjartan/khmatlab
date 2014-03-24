
%% PANEL DEMO 7
%%
%% (a) Create an empty panel object
%% (b) Packing two panels into one figure
%% (c) Plotting into those panels, setting axis labels
%% (d) Making room for titles



%% (a)

% clear or create figure and create root panel covering
% whole figure
clf
p = panel;
p.fontsize = 12;



%% (b)

% form an upper panel A and a lower panel B where A takes 30% of the
% available space, leaving 70% for B
pA = p.pack(30);
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

% set a title on the lower panel
pB.title = ['my panel title'];

% no space is allocated for a title by default in a panel
% layout, so this title bangs into the above xlabel. we can
% fix this explicitly by setting the margin for the panel to
% inclue a space for the title...
pB.axismargin = [15 15 0 5];


