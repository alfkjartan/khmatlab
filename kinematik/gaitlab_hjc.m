function hjc = gaitlab_hjc(mdata, pelvismarkers)
%  hjc = gaitlab_hjc(mdata, pelvismarkers)
% Returns the hip joint center in the same coordinate system as
% mdata using the method from gaitlab (Kit Vaughan). Pelvismarkers
% must be a cell array of strings with the order: ASIS_L, ASIS_R,
% SACRUM.
% hjc = [hjc_l  hjc_r]
 
% Kjartan Halvorsen
% 2004-09-29
  
  
pmarkers = extractmeanmarkers(mdata, pelvismarkers);
asis_l = pmarkers(:,1);
asis_r = pmarkers(:,2);
sacrum = pmarkers(:,3);

% Distance between asis markers
v = asis_l - asis_r;
asis_dist = norm(v);

% Unit vectors in local coordinate system of the pelvis
asis_mid = mean(cat(2,asis_l, asis_r),2);
u = asis_mid - sacrum;
u = u/norm(u);

v = v/norm(v);

w = cross(u,v);

hjc_l = sacrum + 0.598*asis_dist*u + 0.344*asis_dist*v - 0.29*asis_dist*w;
hjc_r = sacrum + 0.598*asis_dist*u - 0.344*asis_dist*v - 0.29*asis_dist*w;

hjc = cat(2, hjc_l, hjc_r);



  