function [gsb, blade] = sensor_frame(md, bfmarkers, fr)
%%  [gsb,blade] = sensor_frame(md, bfmarkers, fr)
%% Computes the transformation taking points in body frame to spatial
%% frame at given frame fr.

%% Kjartan Halvorsen
%% 2012-05-10

mdfr = md;
mdfr{2} = md{2}(fr,:);

origin = extractmarkers(mdfr, bfmarkers.origin);
origin = reshape(origin, 3, length(origin)/3);
origin = mean(origin, 2);

posx = (extractmarkers(mdfr, bfmarkers.posX))';
negx = (extractmarkers(mdfr, bfmarkers.negX))';
posY = (extractmarkers(mdfr, bfmarkers.posY))';
negY = (extractmarkers(mdfr, bfmarkers.negY))';
bladeS = (extractmarkers(mdfr, bfmarkers.blade))';

ex = posx - negx;
ex = ex / norm(ex);
ey = posY - negY;
ey = ey - (ey'*ex)*ex;
ey = ey / norm(ey);
ez = cross(ex,ey);

gsb = eye(4);
gsb(1:3,1) = ex;
gsb(1:3,2) = ey;
gsb(1:3,3) = ez;
gsb(1:3,4) = origin;

blade = ginv(gsb)*cat(1, bladeS, 1);