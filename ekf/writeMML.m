function s=writeMML(km,states,pth,filename,namebase,samplefreq,unit, ...
		    vertical,mfname)
%
% s=writeMML(km,pth,filename,namebase,samplefreq,unit,vertical,mfname)
%
% Writes the kinmodel to the specified file. The file will be in
% the mml (xml) format specified in
% http://www.syscon.uu.se/~kha/modelmotion.dtd 
%
% Input
%    km            ->  kinmodel struct. Contains
%                      km.twists
%                      km.p0
%                      km.gcnames
%                      km.segm_names
%                      km.segm_ends
%                      km.segm_shape
%    states        ->  (dofs x N) matrix
%    pth           ->  path to files to be written
%    filename      ->  name of the mml file to be written
%    namebase      ->  base of the name. Used when writing joint
%                      files
%    samplefreq    ->  Sampling frequency
%    unit          ->  string with length unit (e.g. 'mm')
%    vertical      ->  unit length vector indicating positive
%                      vertical direction (upwards).
%    mfname        ->  name of file to write the marker data

% Based on the writeXML method in the kinmodel class

% Kjartan Halvorsen
% 2003-02-20
%
% Revisions
% 2993-08-25   Changed to more general code (not specific for the
%              eqdist model)

newline=sprintf('\n');

%preamble=['<!DOCTYPE modelmotion SYSTEM "http://www.syscon.uu.se/~kha/modelmotion.dtd">',newline];
preamble='';

samplestr=sprintf('%f',samplefreq);
verticalstr=sprintf('%f ',vertical);

mmstart=['<modelmotion samplefreq="',samplestr,'" unit="',unit,...
  '" verticalup="',verticalstr,'" >',newline];
mmend=['</modelmotion>',newline];

mstart=['<model created="',date, ...
	'" title="equine dist limb">',newline];
mend=['</model>',newline];
mmend=['</modelmotion>',newline];

% The initial marker positions
[mnames, p0] = prepare_mdata(km.p0);

%segm=toXML(km.segm_ends, km.segm_shape, km.twists, states, km.jcs,...
%	   km.gcnames, km.segm_names, fullfile(pth,namebase)) ;
segm=toXML(km.segm_ends, km.segm_shape, km.twists, states, p0, km.gcnames,...
	   km.segm_names, fullfile(pth,namebase)) ;

% Write marker data element
if (nargin==9)
  mdata=['<markerdata radius="16" src="file:',mfname,...
	 '" />',newline];
else
  mdata='';
end
   
fname=fullfile(pth,filename);
fid=fopen(fname,'w');
s=fprintf(fid,'%s%s%s%s%s%s',preamble,mmstart,mstart,segm,mend,mdata,mmend);
fclose(fid);


function [segm, dof, sgms]=toXML(segm_ends, segm_shape, ...
				 twists,states,p0,gcnames,...
				 segmname,fnamebase)
  
  newline=sprintf('\n');

  ssh=cylinderXML(segm_ends{1}, segm_shape{1}, p0{1});

  dof=length(twists{1});

  sgms=1;
  
  segmname{1}
  if (dof>0)
    sjoint=jointXML(twists{1},states(1:dof,:),...
		 gcnames(1:dof),[fnamebase,segmname{1}]);
  else
    sjoint='';
  end
  
  sbr='';
  if ~(length(twists)==1) % Not the last segment
     % Recursive call
     sbr = '';
     for br=2:length(twists)
       [s_br,dof_br,sgms_br] = ...
           toXML(segm_ends{br}, segm_shape{br}, ...
		 twists{br},states(dof+1:end,:),p0{br},...
		 gcnames(dof+1:end),segmname(sgms+1:end),fnamebase);
       sgms = sgms+sgms_br;
       dof = dof+dof_br;
       sbr = [sbr, s_br];
     end
  end
  
  starttag=['<segment title="',segmname{1},'">',newline];
  endtag=['</segment>',newline];

  segm=[starttag,ssh,sjoint,sbr,endtag];


function jstr=jointXML(twists,states,gcnames,fname)

  newline=sprintf('\n');
  starttag=['<jointmodel dof="6">',newline];
  sstate=['<jointstate src="','file:',fname,'"/>',newline];
  endtag=['</jointmodel>',newline];

  jstr=[starttag,sstate,endtag];

  % Write the states to file
  st=rbt2mml12form(twists,states);
  st=st';

  fid=fopen(fname,'w');
  fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',st);
  fclose(fid);
 

  
  
function sh=cylinderXML(segm_ends, segm_shape, p0)
  % In the future change from cylinder to truncated cone
    
  newline=sprintf('\n');

  jc_dist = segm_ends(:,2);
  jc_prox = segm_ends(:,1);
  
  length = norm(jc_dist-jc_prox);
  if (length < 1e-8) 
    sh='';
    return 
  end
  
  radius = mean(segm_shape)*length;
  
  % First the rigid body transformation part.
  e_ax=jc_dist-jc_prox;
  e_ax=e_ax/length;
  rotm=vects2rotation([0;1;0],e_ax);
  rot=rbt2rotation(rotm);
  transl=mean(segm_ends,2);
  rbtstart='<rbt ';
  rbtend=[' />',newline];

  rotstr=sprintf(' %f',rot);
  translstr=sprintf(' %f',transl);

  srbt=[rbtstart, 'rotation="',rotstr, '" translation="', translstr ...
	, '" ', rbtend];

  % Then the geometries
  geomstart=['<shapegeometry>',newline];
  geomend=['</shapegeometry>',newline];

  starttag = '<cylinder ';
  endtag = '/>';

  height=['height="', num2str(length), '" '];
  rad=['radius="', num2str(radius), '" '];

  sr=[starttag,height,rad,endtag,newline];

  rbt=[rotm transl; 0 0 0 1];

  if (~isempty(p0))
    p0_local=inv(rbt)*cat(1,p0,ones(4-size(p0,1),size(p0,2)));
  
    % Write the markers as small spheres
    sspheres='';
    for i=1:size(p0_local,2)
      sspheres=[sspheres,sphereXML(p0_local(1:3,i))];
    end
    smarker = markerXML(p0_local);
  else
    sspheres = '';
    smarker = '';
  end
  
  sgeom=[geomstart,sr,sspheres,geomend];

  
    
  starttag=['<shape>',newline];
  endtag=['</shape>',newline];

  sh=[starttag,srbt,sgeom,smarker,endtag];


  newline=sprintf('\n');

  
  
function s=hoofXML(rbt,p0)
  
  newline=sprintf('\n');

  % First the rigid body transformation part.
  rot=rbt2rotation(rbt);
  transl=rbt(1:3,4)';
  
  rbtstart='<rbt ';
  rbtend=[' />',newline];

  rotstr=sprintf(' %f',rot);
  translstr=sprintf(' %f',transl);

  srbt=[rbtstart, 'rotation="',rotstr, '" translation="', translstr ...
	, '" ', rbtend];

  % Then the geometry
  geomstart=['<shapegeometry>',newline];
  geomend=['</shapegeometry>',newline];

  % The polygonset
  scale=35;
  np=20; % The number of vertices
  [co,coi]=hoofPolygons(scale,np);
  
  polystart='<polygonset ';
  polyend=['/>',newline];

  coord='coord="';
  coord=[coord,sprintf('%f %f %f,',co')];
  coord(end)=[]; % Remove trailing comma
  coord=[coord,'" '];

  coordI='coordInd="';
  coordI=[coordI,sprintf('%f, ',coi')];
  coordI(end)=[];
  coordI=[coordI,'" '];

  p0_local=inv(rbt)*cat(1,p0,ones(4-size(p0,1),size(p0,2)));
  
  % Write the markers as small spheres
  sspheres='';
  for i=1:size(p0_local,2)
     sspheres=[sspheres,sphereXML(p0_local(1:3,i))];
  end
  
  smarker = markerXML(p0_local);
  
  sgeom=[geomstart,polystart,coord,coordI,polyend,sspheres,geomend];

  % Putting it together
  starttag=['<shape>',newline];
  endtag=['</shape>',newline];

  s=[starttag,srbt,sgeom,smarker,endtag];



function smarker = markerXML(p0)
  % the marker positions

  newline=sprintf('\n');
  smarker='';
  sspheres='';

  [m,n]=size(p0);
  mstart='<marker pos="';
  mend=' "/>';
  for i=1:n
    mpos=sprintf(' %f',p0(1:3,i)');
    smarker=[smarker,mstart,mpos,mend,newline];
  end


function xml=sphereXML(p)
  newline=sprintf('\n');

  starttag = '<sphere ';
  endtag = '/>';

  xs=['offsetx="', num2str(p(1)), '" '];
  ys=['offsety="', num2str(p(2)), '" '];
  zs=['offsetz="', num2str(p(3)), '" '];
%  rad='radius="16"';
  rad='radius="20"';

  xml=[starttag,xs,ys,zs,rad,endtag,newline];

  % Putting it together






function  mlf=rbt2mml12form(twists,states)

  N=size(states,2);
  mlf=zeros(N,12);
   
  for t=1:N
    g=eye(4);
    for st=1:length(twists)
      gs=expm(twists{st}*states(st,t));
      g=g*gs;
    end
    mlf(t,:)=reshape(g(1:3,:),1,12);
  end
 
 
  
function [coord,coordInd]=hoofPolygons(scale,np,fname1,fname2)
% Returns a set of coordinates and an indices that defines a set of 
% polygons for modelling the hoof. If the filenames are given, the results
% are written to the respective files.

% Kjartan Halvorsen
% 2000-05-16

  uc=zeros(np+2,3);
  lc=zeros(np+2,3);

  uc(1,:)=[-1 0 1];
  uc(np+2,:)=[-1 0 -1];
  lc(1,:)=[-1 1.8 1.2];
  lc(np+2,:)=[-1 1.8 -1.2];
  ang=pi/2:-pi/(np-1):-pi/2;
  for k=1:length(ang)
    uc(k+1,:)=[0.5+cos(ang(k)) 0 sin(ang(k))];
    lc(k+1,:)=[0.5+1.2*cos(ang(k)) 1.8 1.2*sin(ang(k))];
  end

  coord=[uc;lc]*scale;

  % The indices
  coordInd=(0:np+1)'; % The top polygon
  coordInd=[coordInd;-1;(np+2:2*(np+2)-1)']; % The bottom face
  coordInd=[coordInd;-1;0;np+2;2*(np+2)-1;np+1]; % The back
  for k=0:np
    coordInd=[coordInd;-1;k;k+1;k+np+3;k+np+2];
  end

  if (nargin==4)
    fid=fopen(fname1,'w');
    fprintf(fid,'%f\t%f\t%f\n',coord');
    fclose(fid);
    fid=fopen(fname2,'w');
    fprintf(fid,'%d\n',coordInd);
    fclose(fid);
  end


