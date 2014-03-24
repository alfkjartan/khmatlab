function [plpoints, plnvectors, plnames, cellarea] = get_pliance_model()
  %%  pmod = get_pliance_model()
  %% Returns a model of the pliance handle. The model contains the 
  %% local position and  normal vector of each cell,  organized in
  %% the same order as in the ascii data exported from Novel.
  
  %% Kjartan Halvorsen
  %% 2012-02-09

  circum = 4*14.13*1e-3;
  %%radius = circum/2/pi;
  %%circum = 2*pi*radius;
  radius = 10*1e-3; % 10mm
  length = 226*1e-3; % 226mm
  celllength = length/16;
  %%cellwidth = circum/4;

  cellwidth = 14.13*1e-3;
  cellarea = celllength*cellwidth;
  %cellarea = celllength*cellwidth;
  
  %%keyboard

  %% Let the local coordinate system have z-axis along the long axis of the
  %% cylinder, origin at the center at the start of the mat (row 1), 
  %% and the x-axis pointing towards the middle of cell (1,1). 
  %% That is, cell (1,1) will have position (radius, 0, 0) and 
  %% normal vector (1, 0, 0).
  %% The order of the points are row major, that is, the first 4 points covers the first
  %% circumference.

  plpoints = zeros(3, 16*4);
  plnvectors = zeros(3, 16*4);
  plnames = cell(16*4, 1);

  k=1;
  for row = 1:16
    z = - celllength/2 + row*celllength;
    for col = 1:4
      %%theta = (col-1) * pi/2;
      theta = -(col-1) * cellwidth/radius;
      [x,y,z] = pol2cart(theta, radius, z);
      plpoints(:,k) = [x;y;z];
      [xn,yn,zn] = pol2cart(theta, 1, 0);
      plnvectors(:,k) = [xn;yn;0];%[xn;yn;zn];
      plnvectors(:,k) = plnvectors(:,k)/norm(plnvectors(:,k));
      plnames{k} = sprintf('pos%d%d', row, col);
      k = k+1;
    end
  end

