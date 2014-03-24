function [c, segms, totmass, CoM] = balance2dmodel(varargin)
%  [c, segms] = balance2dmodel(ang1, ang2, ..., mod2d, jointc)
% Returns the squared horizontal distance from the total center of
% mass of the 2d model to the distal end of the distal most
% segment. 
%
% Input
%    angx       ->  The joint angles, corresponding to similar
%                   elements in mod2d and jointc
%    mod2d      ->  Cell array with 2d models. Distal segments
%                   first.
%    jointc     ->  array with 1 or -1 indicating that the rotation axis
%                   comes out of or goes into the plane of the 2d
%                   movement.
% Output
%    segms      ->  Cell array, where each cell
%                   contains a segment struct:
%                   mass   ->  he mass of the segment
%                   CoM    ->  Current position of CoM
%                   prox   ->  Current position of proximal joint
%                   dist   ->  Current position of distal joint

% Kjartan Halvorsen
% 2004-10-05
  
if nargin>0
  mod2d = varargin{end-1};
  jointc = varargin{end};

  segms = cell(size(mod2d));
  
  R = eye(3);
  prox = [0;0;1];

  totmass = 0;
  CoM = [0;0;1];
  
  % assume segments organized as a vertical chain. Distal tip of
  % distal segment at origin.

  for s=1:length(mod2d)
    segm.mass = mod2d{s}.mass;
    segm.dist = prox;
    
    cs = cos(varargin{s});
    sn = sin(varargin{s});
    sgn = sign(jointc(s));
    Rs = [cs -sgn*sn 0; sgn*sn cs 0; 0 0 1];
    R=R*Rs;
  
    segm.CoM = segm.dist + R*[0; mod2d{s}.length - mod2d{s}.CoM(1);0];
    prox = segm.dist + R*[0; mod2d{s}.length;0];
    segm.prox = prox;
    segm.p0 = cat(2, segm.prox', segm.dist');
    segms{s} = segm;
    
    CoM = (totmass*CoM + segm.mass*segm.CoM) / ...
	  (totmass + segm.mass);
    totmass = totmass + segm.mass;
  end

  % The quadratic horizontal distance from CoM to origin (distal
  % end of the chain)
  c = CoM(1)^2;
  
else
  disp('Unit test of balance2dmodel')
  
  tolr = 1e-12;
  % Make some segments
  segm1.mass = 2;
  segm1.length = 3;
  segm1.CoM = [1.5; 0; 1];
    
  segm2.mass = 3;
  segm2.length = 4;
  segm2.CoM = [2.5; 0; 1];
  
  jointc = [1 -1];
  
  ang = pi/3;
  [c,sgms, totmass, CoM] = balance2dmodel(ang, 2*ang, {segm1,segm1}, jointc)
  
  if (abs(c - (segm1.CoM(1)*sin(ang))^2) > tolr | ...
      abs(CoM(2) - segm1.length*cos(ang)) > tolr)
    disp('Unit test failed')
  else
    disp('Unit test passed')
  end    
  
end
