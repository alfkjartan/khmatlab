function nmd=clipmdata(md, tittel)
% This function lets the user choose to clip the data series. 
% Basically, the user chooses two time frames by clicking in the
% plot. After the second click, the data between the two time frames is
% extracted, the rest is discarded, and the data is plotted. 
% The user may undo by hitting 'Escape'.

% Kjartan Halvorsen
% 2001-09-24

if (nargin < 2)
  tittel = 'Select a range by clicking twice in the plot';
end

if (size(md,1)==1) 
   nmd=md;
   return; 
end

hrange=figure;
plot(md);

title(tittel);

[x,y]=ginput(2);

if (length(x==2))
   nmd=md(floor(min(x)):ceil(max(x)),:);
else
   nmd=md;
end

close(hrange)

