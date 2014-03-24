function create6doffile
% create6doffile
% Creates a file with 6dof data and possibly force data from a c3d
% file

% Kjartan Halvorsen
% 2006-02-21

[xmlfile, datahome] = uigetfile('*.xml', 'Open the xml file');

if ~xmlfile
  return
end

cd(datahome)

% Parse the xml file
try
  dom = xmlread(fullfile(datahome, xmlfile));
catch
  error('Could not parse the xml file');
end

% List of nodes
xmlnodes = getChildNodes(dom);

sixdofnodes=getChildNodes(item(xmlnodes,0));

% Find and parse model first
segments = {};
for s = 0:(getLength(sixdofnodes)-1) % Indexed from 0
 subnode = item(sixdofnodes, s);
 if (getNodeType(subnode) == subnode.ELEMENT_NODE)
   nodename = char(getNodeName(subnode));
   if strcmp(nodename, 'model')
     modelnodes = getChildNodes(subnode);
     for sgm = 0:(getLength(modelnodes)-1) % Indexed from 0
       snode = item(modelnodes,sgm);
       if (getNodeType(snode) == snode.ELEMENT_NODE)
	 snodename = char(getNodeName(snode));
	 if strcmp(snodename, 'segment')
	   sgmname = char(getAttribute(snode, 'name'));
	   sgmmarkers = split(char(getAttribute(snode, 'markers')), ...
			      ',');
	   segments = putvalue(segments, sgmname, sgmmarkers);
	 end
       end
     end
   end
 end
end

segments

% Now process the files
for s = 0:(getLength(sixdofnodes)-1) % Indexed from 0
  subnode = item(sixdofnodes, s);
  if (getNodeType(subnode) == subnode.ELEMENT_NODE)
    nodename = char(getNodeName(subnode));
    if strcmp(nodename, 'file')
      file_in = char(getAttribute(subnode, 'name'));
      file_out = char(getAttribute(subnode, 'out'));
      file_extra = char(getAttribute(subnode, 'analogdata'));
      avratio = ...
	  str2double(char(getAttribute(subnode, 'analogvideoratio')));
      disp(['Reading file: ', file_in])
      md = openmocapfile('', fullfile(datahome,file_in));
      if ~isempty(file_extra)
	disp(['Reading file: ', file_extra])
	mdextra = dlmread(fullfile(datahome,file_extra),'\t', 4, 0);
	mdextra = downsample(mdextra, avratio);
      end
      
      sixdofdata = zeros(size(md{2},1), 16*size(segments,1));
      for i = 1:size(segments,1)
	disp(['Computing motion of segment: ', segments{i,1}])
	sixdofdata(:,(i-1)*16+1:i*16) = getMotion(md, segments{i,2});
      end
      
      dlmwrite(fullfile(datahome, file_out), ...
	       cat(2, sixdofdata, mdextra), ...
	       'delimiter','\t');
    end
  end
end

      
