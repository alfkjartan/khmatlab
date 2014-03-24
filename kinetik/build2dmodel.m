function m2d = build2dmodel(proxjc, distjc, antrop_data, ...
			    trackmarkers, markernames)
%  m2d = build2dmodel(prox_jc, dist_jc, antrop_data, ...
%			    trackmarkers, markernames)
%
% Builds a single 2d segment model for inverse dynamics.
% 
% Input
%    prox_jc          ->  Proximal joint center (3 x 1) 3d point
%    dist_jc          ->  Distal joint center (3 x 1) 3d point.
%    antrop_data      ->  Struct containing the antropometric data
%                         mass  -> The mass of the segment
%                         CoM   -> The position of center of mass
%                                  given in procent of the distance
%                                  from the proximal to the distal
%                                  joint center. If empty, the the
%                                  last two fields must be given.
%                         I     -> Moment of inertia. If empty,
%                                  then the next two fields must be
%                                  provided.
%                         length  -> Optional. Otherwise, the
%                                    distance between prox_jc and
%                                    dist_jc is used.
%                         prox_circumf -> Circumference at the
%                                         proximal end.
%                         dist_circumf -> Circumference at the
%                                         distal end.
%    trackmarkers  ->  Position of markers used for tracking. 
%                        (3 x 1) matrix.
%    markernames     ->  Cell array with names of tracking markers.
% Output
%    m2d             <-  Struct representing the 2D segment model.
%                        mass   ->  The mass of the segment  
%                        length ->  The length of the segment  
%                        CoM    ->  The center of mass of the segment  
%                        I      ->  The moment of inertia of the segment  
%                        p0     ->  Position of markers on segm in local
%                                   coordinate system. (1 x 3m) row
%                                   vector. The positions must lie
%                                   in the z=1 plane. This is a
%                                   trick to be able to use the
%                                   points as homogenous points in 2D 
%                        markernames  ->  Cell array with names
%                                         of tracking markers. 
%                        

% Kjartan Halvorsen
% 2004-03-12


% -----------------------------------------------------------------
% The average density of the body. Hardcoded constant
bodydensity = 1.06e3; % kg/m^3 
% -----------------------------------------------------------------

% Define transformation which takes points in the local coordinate
% system to the global.

d_gl = proxjc(1:2);
e_x = distjc-proxjc;
e_x = e_x / norm(e_x);
e_x_local = [1;0;0];

G_gl = soder2d(cat(1,[0 0 1 0], [d_gl' (d_gl+e_x(1:2))']));
G_lg = soder2d(cat(1, [d_gl' (d_gl+e_x(1:2))'], [0 0 1 0]));

if 0
  tmp = cross(e_x_local, e_x);
  th = sign(tmp(3))*acos(e_x'*e_x_local);


  R_gl = [cos(th) -sin(th); sin(th) cos(th)];
  G_gl = cat(1, cat(2, R_gl, R_gl'*d_gl),...
	     [0 0 1]);

  % The inverse transformation
  G_lg = cat(1, cat(2, R_gl', -d_gl),...
	     [0 0 1]);
end



% Convert the trackmarkers to local coordinates.
p0=G_lg * trackmarkers;


m2d.p0 = (p0(:))';

if isfield(antrop_data,'length')
  m2d.length = antrop_data.length;
else
  m2d.length = norm(distjc(1:2)-proxjc(1:2));
end

% Calculate moment of inertia if needed. 
if (~isfield(antrop_data,'I') | isempty(antrop_data.I))
  if (~isfield(antrop_data,'prox_circumf') | ...
      isempty(antrop_data.prox_circumf) | ...
      ~isfield(antrop_data,'dist_circumf') | ...
      isempty(antrop_data.dist_circumf))
    error(['The circumferences at proximal and distal end of the segment' ...
	   ' are needed'])
  end

  radprox = antrop_data.prox_circumf/2/pi;
  raddist = antrop_data.dist_circumf/2/pi;
  
  [m, r0, I0, I] = coneinertia(m2d.length, radprox, raddist);
  
  m2d.I = I*bodydensity;
  m2d.mass = m*bodydensity;
  
  m2d.CoM = [r0; 0; 1];
else
  m2d.I = antrop_data.I;
  m2d.mass = antrop_data.mass;
  m2d.CoM = [m2d.length*antrop_data.CoM/100; 0; 1];
end

m2d.markernames = markernames;

m2d.prox = proxjc;
m2d.dist = distjc;


  




  