function [vel, acc, angvel, angacc] = segmentvelacc2d(mdata, ...
						  model2d, ...
						  plane_norm, vertical)
% [vel, acc, angvel, angacc] = segmentvelacc2d(mdata, model2d, ...
%                                              plane_norm, vertical)
%
% Computes the velocity and angular accelleration of a segment in
% 2D.
% 
% Input
%    mdata        ->  mdata struct
%    model2d      ->  2d model struct. Fields used in this function:
%                     CoM    ->  The center of mass of the segment  
%                     p0     ->  Position of markers on segm in local
%                                coordinate system. (1 x 3m) row
%                                vector. The positions must lie
%                                in the z=1 plane. This is a
%                                trick to be able to use the
%                                points as homogenous points in 2D 
%                     markernames  ->  Cell array with names
%                                      of tracking markers. 
%    plane_norm   ->  Unit vector normal to 2D plane.
%    vertical     ->  Unit vector in the vertical direction of the 2D plane.
% Output
%    vel          <-  Velocity of CoM
%    acc          <-  Acceleration of CoM
%    angvel       <-  Angular velocity
%    angacc       <-  Angular acceleration
  
% Kjartan Halvorsen
% 2004-03-16

if (nargin > 0)

  fs = str2double(getvalue(mdata{1},'FREQUENCY'));
  %dt = 1/fs;

  md = extractmarkers(mdata, model2d.markernames);
  
  if (nargin > 2)
    md = project2plane(md, plane_norm, vertical);
  end
  
  [nfr, m3] = size(md);

  pos = zeros(nfr,3);
  ang = zeros(nfr,1);
  angvel = zeros(nfr,1);

  tol = 1e-9; % tolerance for determinant of rotation matrix

  p0soder = model2d.p0;
  p0soder(3:3:end)=[];
  
  mm = m3/3;
  p0 = reshape(model2d.p0,3,mm);
  p0mean = mean(p0,2);
  axx = diff(p0(:,1:2),1,2);
  axx = axx/norm(axx);
  
  for i=1:nfr

    use_soder2d = 0;
    
    if use_soder2d
      ppi = md(i,:);
      ppi(3:3:end) = [];
      G2di = soder2d(cat(1, p0soder, ppi));

      pos(i,:) = (G2di*model2d.CoM)';
    
      % Find the incremental rotation angle
      ppimin1 = md(max(i-1,1),:);
      ppimin1(3:3:end) = [];
      G2dimin1 = soder2d(cat(1, ppimin1,p0soder));
      G = G2dimin1*G2di;
    
      angvel(i) = ...
	  sign(G(2,1))*(acos(G(1,1)) + acos(G(2,2)))/2*fs;

      ang(i) = ...
	  sign(G2di(2,1))*(acos(G2di(1,1)) + acos(G2di(2,2)))/2;
      
    else
      % Find the position of the center of mass
    
      ppi = reshape(md(i,:),3,mm);
      ppimean = mean(ppi,2);
      displ = ppimean - p0mean;
    
      axxi = diff(ppi(:,1:2),1,2);
      axxi = axxi/norm(axxi);
    
      ph = sign([0 0 1]* cross(axx,axxi))*acos(axx'*axxi);
      ang(i) = ph;
      pos(i,:) = (ppimean + [cos(ph) -sin(ph) 0; sin(ph) cos(ph) 0; ...
			0 0 1]*(model2d.CoM-p0mean))';
    
      ppimin1 = reshape(md(max(i-1,1),:),3,mm);
      axximin1 = diff(ppimin1(:,1:2),1,2);
      axximin1 = axximin1/norm(axximin1);

      angvel(i) = ...
	  sign([0 0 1]* cross(axximin1,axxi))*acos(axximin1'*axxi)* ...
	  fs;
    end
  end

  ang;
  

  vel = centraldiff(pos,fs);
  acc = centraldiff(vel,fs);

  vel = vel';
  acc = acc';

  use_increment_angles = 1;
  if ~use_increment_angles
    angvel = centraldiff(ang,fs);
  end

  angacc = centraldiff(angvel,fs);

elseif (nargin==0) % Unit test
  
  % Set up test case
  tlngth = 80;
  pp = zeros(3,3,tlngth);
  phi = (linspace(0,pi,tlngth)).^3;

  phi';
  p0 = cat(1,randn(2,3),ones(1,3));
  
  for i=1:tlngth
    pp(:,:,i) = [cos(phi(i)) -sin(phi(i)) 0
		sin(phi(i)) cos(phi(i)) 0
		0 0 1]*p0;
  end
  
  mnames = {'p1','p2','p3'};
  
  attr = {'FREQUENCY', '1'
	  'MARKER_NAMES',mnames};
  mdata = {attr, (reshape(pp,9,tlngth))'};
  
  
  mod2d.CoM = ones(3,1);
  mod2d.p0 = reshape(p0,1,9);
  mod2d.markernames = mnames;
  
  plane_norm = [0;0;1];
  vertical = [0;1;0];

 % Call the function
 [vvel, acc, angvel, angacc] = segmentvelacc2d(mdata, ...
					      mod2d);
 
 % Test angvel och angacc.
 
 angvel_exp = centraldiff(phi',1);
 angacc_exp = centraldiff(angvel_exp,1);

 figure
 clf
 subplot(211)
 plot(1:tlngth,angvel_exp,'b', 1:tlngth,angvel,'r')
 subplot(212)
 plot(1:tlngth,angacc_exp,'b', 1:tlngth,angacc,'r')
  
  
  
end
