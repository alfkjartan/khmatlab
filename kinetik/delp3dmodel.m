function delp3d = delp3dmodel
% delp3d = delp3dmodel
% 
% Loads (or generates) a 3d model of the lower límbs based on delps
% model (Delp et al. " An interactive graphics-based model of the
% lower extremity to study orthopaedic surgical procedures". IEEE
% tr Biomed Engn, Vol 37, 1990) 
%
% The function returns a struct containing the following fields:
%    segments   <-  A set of segments. These have the fields:
%                   name      <-  string
%                   coordsys  <-  (4 x 4) local coordinate system
%    joints     <-  A set of joints. These have the fields:
%                   name      <-  string
%                   proximal  <-  the proximal segment
%                   distal    <-  the distal segment
%                   tws       <-  a set (cell array) of twists
%    muscles    <-  A set of muscles. These have the fields:
%                   name      <-  string
%                   coordsys  <-  (4 x 4) local coordinate system

  
%                   
% Segments
world_segm = 
  
% Joints

pelvis.proximal = 
% Muscles
  
