function [tws, p0, gcnames, jc, segmnames, CoM, radius, mass, g0, inertia, objectframe, objectcenter]=build_model(varargin)
% Builds and returns a kinematic model (not a matlab object). Based
% on build_km
%
% Usage: 
%  [tws, p0, gcnames, jc, segmnames, CoM, radius, mass] = build_model(sgm1,sgm2,...)
% Input
%   sgmi   ->   Struct containing the information needed to 
%               define a segment. The fields of the struct are
%                 name       -> string
%                 localframe -> (4x4) matrix defining the local 
%                               coordinat system of the segment.
%                               Used for defining the orientation 
%				of the joint axes, as well as the 
%				center of the joint (origo). 
%                 shapeframe -> (4x4) matrix defining the local 
%                               coordinat system for the shape. 
%				The y-direction is the main axis of
%				the segment.
%                 dof        -> (2 x 1) cell array describing the
%			        degrees of freedom  
%				of the ith joint. The first element gives the 
%				dof for rotation. The second gives the dof for
%				translation. By definition, the axis
%				of rotation, as well as the directions of
%				translation are the local x, y, and z axis. 
%				Example: dof={[3 2] [1 2 3]} gives a jointmodel
%		  		with translation in all directions
%				followe by rotation around the z-axis, then
%				the y-axis.
%                 markers    -> a hash containing the name and the
%		                position  of the markers associated
%				with the segment.
%                 valmarkers -> Similar to markers, and contains
%		                additional markers such as bone-fixed markers
%                               used for validating the calculated motion.
%                 mass       -> mass of the segment
%                 moment_of_inertia      -> 3x3 moment of inertia wrt its CoM and in the local
%                               coordinate system defined by g0
%                 g0         -> rigid trf defining the local coordinate system for which
%                               moment_of_inertia is defined. g0 transforms vectors in local
%                               coordinate system to global at reference configuration (th=0).
%                 generalized_inertia  -> 6x6 matrix [m*I, 0; 0; mom_of_inertia]
%                 objectframe  -> 4x4 matrix giving the local frame of the manipulated object.
%                                  Should only be given for end segment.
%
% Output
%    tws         <-   An nested cell array of twists
%    p0          <-   Nested cell array with reference marker positions
%    gcnames     <-   Names of generalized coordinates (joint angles)
%    jc          <-   Nested cell array with joint centra.
%    segmnames   <-   Names of generalized coordinates (joint angles)
%    CoM         <-   Nested cell array with centers of mass 
%                     (if given for the segments).
%    radius      <-   Nested cell array with radii for the shapes
%                     (if given for the segments).
%    mass        <-   Vector containing the mass for each segment.
%    g0          <-   Nested array of rigid trf. 
%    inertia     <-   Nested array of generalized inertia.
%    objectframe   <-   Nested array of object_frames.
%    objectcenter   <-   Nested array of object_frames.
%
%

% ------------------------------------------------------------
% Kjartan Halvorsen
% 2002-12-10
%
%
% Revisions
% 2004-03-22   Added output jc - joint centers. Assumed to be the
%              center of the local coordinate system
% 2009-06-25   Changed output mass to be nested array. 
%              Also changed outputs CoM and jc to have same
%              structure as the marker positions.
% 2013-06-05   Changed to work with g0, a local coordinate system with center at the CoM
%              and inertia, nested array of 6x6 inertia matrices.
% ------------------------------------------------------------

nsegm=nargin; % The number of segments.

% Start with the last, most distal segment

% The local coordinate systems are oriented so that y is in the
% direction of the axis, x is in the lateral-medial direction and z 
% is in the direction normal to the other two (posterior-anterior
% or anterior-posterior).


for s=nsegm:-1:1

  if (s==nsegm)
    tws=cell(1);
    p0=cell(1);
    CoM=cell(1);
    radius=cell(1);
    jc=cell(1);
    gcnames = {};
    mass = {};
    segmnames = {};
    g0 = {};
    inertia = {};
    objectframe = {};
    objectcenter = {};
  else
    twsbr=tws;
    tws=cell(2,1);
    tws{2}=twsbr;
    p0br=p0;
    p0=cell(2,1);
    p0{2}=p0br;
    CoMbr = CoM;
    CoM=cell(2,1);
    CoM{2}=CoMbr;
    radbr = radius;
    radius=cell(2,1);
    radius{2}=radbr;
    jcbr = jc;
    jc=cell(2,1);
    jc{2}=jcbr;
    g0br = g0;
    g0 = cell(2,1);
    g0{2} = g0br;
    inertiabr = inertia;
    inertia = cell(2,1);
    inertia{2} = inertiabr;
    objectframebr = objectframe;
    objectframe = cell(2,1);
    objectframe{2} = objectframebr;
    objectcenterbr = objectcenter;
    objectcenter = cell(2,1);
    objectcenter{2} = objectcenterbr;
  end

  segm=varargin{s};

  if ~isfield(segm,'name')
    segm.name = 'unnamed';
  end
  
  g=segm.localframe;
  center=g(1:3,4);

  rdof=segm.dof{1};
  tdof=segm.dof{2};
  jm=cell(length(rdof)+length(tdof),1);
  k=0;
  for trnsl=tdof;
    k=k+1;
    w=zeros(3,1);
    v=g(1:3,trnsl);
    jm{k}=hat([v;w]);
  end
  for rot=rdof
    k=k+1;
    w=g(1:3,rot);
    v=-cross(w,center);
    jm{k}=hat([v;w]); 
  end

  tws{1}=jm;
  
  if isfield(segm, 'markers')
    p0{1}=segm.markers;
  end

  if isfield(segm, 'CoM')
    CoM{1} = {[segm.name,'_CoM'], segm.CoM};
  end

  if isfield(segm, 'radius')
    radius{1} = {[segm.name,'_radius'], segm.radius};
  end
  
  if isfield(segm, 'mass')
    %mass(s) = segm.mass;
    mass = cat(1, {segm.mass}, mass);
  else
    mass = cat(1, {0}, mass);
  end    
  
  jc{1} = {[segm.name,'_jc'], center};
  
  gcnames=cat(1,segm.states,gcnames);
  segmnames=cat(1,{segm.name},segmnames);

  if isfield(segm, 'g0')
    g0{1} = segm.g0;
    %%disp([segm.name, ' got g0']), g0
  end

  if isfield(segm, 'generalized_inertia')
    inertia{1} = segm.generalized_inertia;
    %%disp([segm.name, ' got inertia']), segm.generalized_inertia
  end

  if isfield(segm, 'object_frame')
    objectframe{1} = segm.object_frame;
    %%g0{1} = segm.object_frame;
    jc{1} = {[segm.name,'_jc'], center};
    objectcenter{1} = {[segm.name, '_objcenter'], segm.object_frame(1:3,4)};
  end
end

