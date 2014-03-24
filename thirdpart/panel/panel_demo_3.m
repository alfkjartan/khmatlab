
%% PANEL DEMO 3
%%
%% (a) Create an empty panel object
%% (b) Packing some panels
%% (c) Plotting into those panels, setting axis labels
%% (d) Adjusting the layout using different units



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

% switch to inches (default is mm)
p.units = 'in';

% observe the margin of pA
disp(['pA.axismargin'])
disp(pA.axismargin)

% set the margin of pA
pA.axismargin = [1 0.59055 0 0];

% switch to mm
p.units = 'mm';

% observe the margin of pA
disp(['pA.axismargin'])
disp(pA.axismargin)
