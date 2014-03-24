function mp = project_on_plane(p, a1, a2)
%%  mp = project_on_plane(p, a1, a2)
%% Projects points p on the plane whos normal is given by the line goint
%% through the points a1 and a2

v = a2-a1;
vnorm = sqrt(sum(v.^2, 2));

