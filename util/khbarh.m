function h=khbarh(dta)
%h=khbarh(dta)
% Similar to barh, except it writes the value of each bar
% on the bar.

% Kjartan Halvorsen
% 2008-11-07

barh(dta);

range = max(dta) - min(dta);

for n = 1:size(dta,1)
  text(dta(n) + range/40, n, sprintf('%0.5f', dta(n)));
  
end