
close all
figure(1)
clf
h = surf (peaks (25));
set (h, "facealpha",0.6)
xlabel ("my xlabel")
ylabel ("my ylabel")
title ("my title")
papersize = [8.5, 11.0];
figuresize = [6, 5];
set (gcf, "papersize", papersize,
          "paperposition", [(papersize-figuresize)/2, figuresize])

print (gcf, "-dtikz", "-color", "gnuplot-tikz-lua.tex")

if true
  more off
  disp ("running latex")
  [status, output] = system ("latex lua_test.tex");
  disp ("running dvipdf")
  [status, output] = system ("dvipdf lua_test.dvi");
  disp ("running preview")
  [status, output] = system ("preview lua_test.pdf");
  more on
end
