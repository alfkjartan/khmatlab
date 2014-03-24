
%% PANEL DEMO 6
%%
%% (a) Create an empty panel object
%% (b) Packing some panels
%% (c) Plotting into those panels, setting axis labels
%% (d) Using engineering scales



%% (a)

% clear or create figure and create root panel covering
% whole figure
clf
p = panel;



%% (b)

% form an upper panel A and a lower panel B where A takes 30% of the
% available space, leaving 70% for B
pA = p.pack(50);
pB = p.pack();



%% (c)

% plot into panel A
select(pA)
x = 1e-6 * (-10:10);
y = (x * 1e7).^2 * 8;
plot(x, y, 'r')

% set axis labels through panel, for control over metrics
pA.xlabel = 'foo';
pA.ylabel = 'bar';

% plot into panel B
select(pB)
plot(x, y, 'r')

% set axis labels through panel, for control over metrics
pB.xlabel = 'foo';
pB.ylabel = 'bar';



%% (d)

% in this case, our axis labels are a bit wacky, and we've
% overflowed the standard margins as a result

% we can fix both things by using engineering scales on both
% axes of the plot (we'll adjust the lower plot only)
pB.xscale = '?';
pB.yscale = '?';

% but now we're left with no indication that the scales are
% not unitary! we can fix this by using the symbol "$" in
% the corresponding labels
pB.xlabel = 'foo ($m)';
pB.ylabel = 'bar ($m/s)';

% additionally, if we wanted, though it's a little
% unconventional, we can have the tick labels suffixed
% appropriately
% pB.yscale = '?$';



