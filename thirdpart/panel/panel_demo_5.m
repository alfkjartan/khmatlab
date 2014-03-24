
%% PANEL DEMO 5
%%
%% (a) Create an empty panel object
%% (b) Packing some panels
%% (c) Plotting into those panels, setting axis labels
%% (d) Adjusting the root margin



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

% panels generally overflow slightly into their parent to
% the top and to the right, since axis number labels stick
% out, just a little bit. we ignore this throughout, since
% it makes layout much easier, but at the root, we don't
% want these labels to go off-figure (also we want a little
% space top and right between the axes and the edges of the
% figure). to see this in action, let's turn *off*
% "rootmargin":
disp('current rootmargin:')
disp(p.rootmargin)
disp('removing rootmargin completely:')
p.rootmargin = [0 0 0 0];
disp(p.rootmargin)
disp('press a key to continue...')
pause

% see, we've lost a bit of those labels and the axes are
% hard up against the top and right edges? so "rootmargin"
% is an extra margin used only by the root panel, and is
% usually used to account for this overflow, thus:
disp('putting back the 5mm top & right rootmargin:');
p.rootmargin = [0 0 5 5]; % 5mm top and right
disp(p.rootmargin)
disp('press a key to continue...')
pause

% however, you can use it as you see fit...
disp('putting in some other rootmargin:');
p.rootmargin = [10 20 30 40];
disp(p.rootmargin)


