
%% PANEL DEMO 4
%%
%% (a) Create an empty panel object
%% (b) Packing multiple panels, using different edges, and
%% using absolute positioning as well as relative
%% (c) Plotting into those panels, setting axis labels
%% (d) Explicit render
%% (e) Export to file



%% (a)

% clear or create figure and create root panel covering
% whole figure
tic
clf
p = panel;

% on complex panels, performance can be improved by turning
% off autorender (this avoids re-rendering the whole figure
% every time any setting is changed). leave autorender on
% whilst you're experimenting, so you'll see the effect of
% any changes you make immediately.
p.autorender = false;



%% (b)

% form an upper panel A and a lower panel B
pA = p.pack(30);
pB = p.pack();

% pack pB on left-hand edge
pB.edge = 'l';

% pack two panels into B
pC = pB.pack(20);
pD = pB.pack();

% pack D with some absolutely-positioned panels
pD1 = pD.pack([0.0 0.3 0.3 0.7]);
pD2 = pD.pack([0.3 0.7 0.7 0.3]);
pD3 = pD.pack([0.7 0.0 0.3 0.7]);
pD4 = pD.pack([0.0 0.0 0.7 0.3]);

% and one in the middle that we can pack more
% relatively-positioned panels into
pD5 = pD.pack([0.3 0.3 0.4 0.4]);

% pack those extra panels into the D central panel
pE = pD5.pack(50);
pF = pD5.pack(50);



%% (c)

% plot into panel
select(pA)
x = -100:100;
plot(x,x.^2/1000,'r')
pA.xlabel = 'x is great';
pA.ylabel = 'y is better';

% plot into panel (no plot, just show axis)
select(pC)

% plot into panel
select(pD1)
plot(randn(1000,3))
pC.xlabel = 'foo';
pC.ylabel = 'bar';

% plot into panel
select(pD2)
plot(randn(1000,3))
pD2.xlabel = 'foo';
pD2.ylabel = 'bar';

% plot into panel
select(pD3)
plot(randn(1000,3))
pD3.xlabel = 'donkey';
pD3.ylabel = 'kong';

% plot into panel
select(pD4)
plot(randn(1000,3))
pD4.xlabel = 'foo';
pD4.ylabel = 'bar';

% plot into panel (no plot, just show axis)
select(pE)
select(pF)



%% (d)

% since we turned off autorender, we have to render the
% layout manually once we've finished building it
render(p);
toc



%% (e)

% export the layout
export(p, '-w200', '-h200', '-rf', '-opng', 'panel_demo.png');

