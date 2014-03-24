function [Fprox, Mprox, Finert, Minert] = invdynsol2d(segm, kinem, forces)
  
%   [Fprox, Mprox] = invdynsol2d(segm, kinem, forces)
%
% Solves the 2d inverse dynamics problem. The reference frame
% considered is fixed in the segment with the origin coinciding
% with the proximal joint center and the x-axis pointing towards
% the distal joint center.
%
% Input
%    segm          ->  struct with the following fields:
%                mass          ->  The mass of the segment  
%                length        ->  The length of the segment  
%                CoM           ->  The center of mass of the segment  
%                I             ->  The moment of inertia of the segment  
%                p0            ->  Position of markers on segm in local
%                                  coordinate system. (1 x 3m) row
%                                  vector. The positions must lie
%                                  in the z=1 plane. This is a
%                                  trick to be able to use the
%                                  points as homogenous points in 2D 
%    kinem         ->  struct with kinematic data:
%                p             ->  Position of points on the
%                                  segment. (1 x 3m) row
%                                  vector. The positions must lie
%                                  in the z=1 plane.
%                com_acc       ->  The accelleration of the center
%                                  of mass. (3 x 1) vector in the
%                                  z=1 plane.
%                ang_acc       ->  The angular accelleration
%    forces        ->  struct with data on external forces: 
%                Fdist         ->  The contact force at the distal
%                                  joint. Must be homogenous vector
%                                  in 2D. 
%                Mdist         ->  The moment at the distal joint
%                g             ->  The direction of gravity
%                Fexternal     ->  Optional. Any other external forces. 
%                                  (3 x n) matrix where each column contains
%                                  the homogenous force vector.
%                pexternal     ->  Optional. Mandatory if Fexternal given.
%                                  (3 x n) matrix where each column contains
%                                  the application point (in
%                                  homogenous coordinates) of the force vector.
%                Mexternal     ->  Optional. Any free external
%                                  moments acting on the segment.
% Output
%    Fprox         <-  The contact force at the proximal joint center
%    Mprox         <-  The moment at the proximal joint center
%

% Kjartan Halvorsen
% 2004-03-10

if (nargin == 3)

  % Compute the rotation matrix and displacement.

  tol = 1e-9; % tolerance for determinant of rotation matrix

  % Find the position of the center of mass
  p0 = segm.p0;
  p0(3:3:end)=[];
  pp = kinem.p;
  pp(3:3:end)=[];
  
  G = soder2d(cat(1, pp, p0)); % Takes global coordinates to local
  Ginv = soder2d(cat(1,p0,pp));

  % Transform center of mass accelleration and
  % external forces
  com_acc = G*kinem.com_acc;

  if ( isfield(forces,'Fexternal') & ~isempty(forces.Fexternal) )
    Fexternal = G*cat(2,segm.mass*9.82*forces.g, ...
		      forces.Fdist,forces.Fexternal);
    pexternal = cat(2, segm.CoM, [segm.length;0;1], ...
		    G*forces.pexternal);
  else
    Fexternal = G*cat(2,segm.mass*9.82*forces.g, forces.Fdist);
    pexternal = cat(2,segm.CoM, [segm.length;0;1]);
  end

  if ( isfield(forces,'Mexternal') & ~isempty(forces.Mexternal) )
    forcesMexternal = forces.Mexternal;
  else
    forcesMexternal = 0;
  end
  
  % The equilibrium equations:

  % Force
  Finert = Ginv*segm.mass*com_acc;
  
  Fprox = Finert - Ginv * sum(Fexternal,2);

  % Moment

  % Externprojekt\bromsman\squat_pilot\al moments
  Mexternal = pexternal(1,:).*Fexternal(2,:) ...
      - pexternal(2,:).*Fexternal(1,:);

%  keyboard
  
%mdist=forces.Mdist

  if ( length(segm.I(:)) > 1 ) % 3 x 3 matrix
    I = segm.I(3,3);
  else
    I = segm.I;
  end

  Minert = I*kinem.ang_acc;
  
  Mprox = Minert - forces.Mdist ...
	  - sum(Mexternal) - sum(forcesMexternal);

%  keyboard
  
elseif (nargin == 0) % Unit test. 
  % Generate a test case. Run function. Check results.

%  segm1.mass = 2;
  segm1.mass = 0;
  segm1.length = 3;
  segm1.CoM = [1.5; 0; 1];
  segm1.I = 0.5*segm1.mass*segm1.length^3;
  segm1.p0 = [0 0 1 3 0 1];
  
  kinem1.p = [1 3 1 1 0 1];
  kinem1.com_acc = [0 0 0]';
%  kinem1.ang_acc = pi;
  kinem1.ang_acc = 0;
  kinem1.ang_acc = 0;
  
  forces1.Fdist = [100; 0; 0];
%  forces1.Fdist = [0; 0; 0];
  forces1.Mdist = 100;
  forces1.g = [0; -1; 0];
  forces1.Fexternal = [0; 20; 0];
  forces1.pexternal = [1; 0; 1];
  forces1.Mexternal = [200];
  
  segm2.mass = 0;
%  segm2.mass = 2;
  segm2.length = 4;
  segm2.CoM = [2.5; 0; 1];
  segm2.I = 0.5*segm2.mass*segm2.length^3;
  segm2.p0 = [0 0 1 4 0 1];
  
  kinem2.p = [-3 3 1 1 3 1];
  kinem2.com_acc = [0 0 0]';
  kinem2.ang_acc = pi/2;
  kinem2.ang_acc = 0;
  
  [F1prox, M1prox] = invdynsol2d(segm1,kinem1,forces1);
  forces2.Fdist = -F1prox;
  forces2.Mdist = -M1prox; 
  forces2.g = [0; -1; 0];
  forces2.Fexternal = [0; 20; 0];
  forces2.pexternal = [1; 0; 1];
  forces2.Mexternal = [];
  
  [F2prox, M2prox] = invdynsol2d(segm2,kinem2,forces2);
  
  
  % Same computations, but reversed chain
  segm2b = segm2;
  segm2b.CoM = [1.5; -0; 1];
  
  kinem2b = kinem2;
  kinem2b.p = [1 3 1 -3 3 1];

  forces2b = forces2;
  forces2b.Fdist = F2prox;
  forces2b.Mdist = M2prox; 

  segm1b = segm1;
  segm1b.CoM = [1.5; -0; 1];
  
  kinem1b = kinem1;
  kinem1b.p = [1 0 1 1 3 1];

  [F2bprox, M2bprox] = invdynsol2d(segm2b,kinem2b,forces2b);

  forces1b = forces1;
  forces1b.Fdist = -F2bprox;
  forces1b.Mdist = -M2bprox;

  forces1a = forces1b;
  forces1a.Fdist = F1prox;
  forces1a.Mdist = M1prox;
  
  [F1bprox, M1bprox] = invdynsol2d(segm1b,kinem1b,forces1b);
  [F1aprox, M1aprox] = invdynsol2d(segm1b,kinem1b,forces1a);

  disp('Unit test of invdynsol2d')

  t1ok = (forces1.Fdist == F1bprox);
  if t1ok
    disp('Test1:  OK');
  else
    disp(['Test1: Failed. Error in distal force of segment 1'])
  end
  
  t2ok = (forces1.Mdist == M1bprox);
  if t2ok
    disp('Test2:  OK');
  else
    disp(['Test2: Failed. Error in distal moment of segment 1'])
    disp(['       Expected ', num2str(forces1.Mdist), ...
	  '    found ', num2str(M1bprox)])
  end
  
  t3ok = (forces2.Fdist == F2bprox);
  if t3ok
    disp('Test3:  OK');
  else
    disp(['Test3: Failed. Error in distal force of segment 2'])
  end
  
  t4ok = (forces2.Mdist == M2bprox);
  if t4ok
    disp('Test4:  OK');
  else
    disp(['Test4: Failed. Error in distal moment of segment 2'])
    disp(['       Expected ', num2str(forces2.Mdist), ...
	  '    found ', num2str(M2bprox)])
  end

    t1ok = (forces1.Fdist == F1aprox);
  if t1ok
    disp('Test5:  OK');
  else
    disp(['Test5: Failed. Error in distal force of segment 1'])
  end
  
  t2ok = (forces1.Mdist == M1aprox);
  if t2ok
    disp('Test6:  OK');
  else
    disp(['Test6: Failed. Error in distal moment of segment 1'])
    disp(['       Expected ', num2str(forces1.Mdist), ...
	  '    found ', num2str(M1bprox)])
  end

%  keyboard

end

