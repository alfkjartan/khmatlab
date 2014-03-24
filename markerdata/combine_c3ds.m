function combine_c3ds(varargin)
% Usage:
% combine_c3ds(c3d_out_file, c3d_in_file_1, c3d_in_file_2,...)
%
% Combines the analog data (so far) in the different c3d files, and
% writes to a new c3d file. The data series can be of different
% length, but it is assumed that they start synchronised.
%

% Kjartan Halvorsen
% 2005-01-20

file_out = varargin{1};

analogs = {};
parameters = {};

channels = 0;
nfrs = 0;
for fil = 2:length(varargin)
  [mrks, vidFrRt, analogs1, analogFrRt, events, parameterGr,...
  caminfo,residerr] = readC3D(varargin{fil});
  
  analogs{fil-1} = analogs1;
  parameters{fil-1} = parameterGr;
  
  channels = channels + size(analogs1,2);
  if size(analogs1,1) > nfrs
    nfrs = size(analogs1,1);
  end
end

parameterGr = parameters{1};

[slsk, angr] = intersect([parameterGr.name], 'ANALOG');
[slsk,usedind] = intersect([parameterGr(angr).Parameter.name],...
			   'USED');
[slsk,lblind] = intersect([parameterGr(angr).Parameter.name],...
			   'LABELS');
[slsk,dind] = intersect([parameterGr(angr).Parameter.name],...
			   'DESCRIPTIONS');
[slsk,scind] = intersect([parameterGr(angr).Parameter.name],...
			   'SCALE');
[slsk,ofind] = intersect([parameterGr(angr).Parameter.name],...
			   'OFFSET');

[slsk, fpgr] = intersect([parameterGr.name], 'FORCE_PLATFORM');
if fpgr
  [slsk,fpuseind] = intersect([parameterGr(fpgr).Parameter.name],...
			      'USED');
  parameterGr(fpgr).Parameter(fpuseind).data = length(analogs);
  
  [slsk,cornind] = intersect([parameterGr(fpgr).Parameter.name],...
			     'CORNERS');

  [slsk,orind] = intersect([parameterGr(fpgr).Parameter.name],...
			     'ORIGIN');

  [slsk,fpchind] = intersect([parameterGr(fpgr).Parameter.name],...
			     'CHANNEL');

  [slsk,zrind] = intersect([parameterGr(fpgr).Parameter.name],...
			     'ZERO');

  [slsk,tpind] = intersect([parameterGr(fpgr).Parameter.name],...
			     'TYPE');

end

% Now combine the information
analog = zeros(nfrs,channels);
chind = 1;
for fil = 1:length(analogs)
  analog(1:size(analogs{fil},1), chind:chind+size(analogs{fil},2)-1) ...
      = analogs{fil};
  
  chind = chind + size(analogs{fil},2);
  if fil>1 
    % Fix the parameters
    parameterGr(angr).Parameter(usedind).data = chind-1;

    lbls = parameters{fil}(angr).Parameter(lblind).data;
    for ch = 1:size(analogs{fil},2)
      lbls{ch}(3) = int2str(fil);
    end
    parameters{fil}(angr).Parameter(lblind).data = lbls;
    parameters{fil}(angr).Parameter(dind).data = lbls;
    
    parameterGr(angr).Parameter(lblind).data = ...
	cat(2,parameterGr(angr).Parameter(lblind).data,...
	    parameters{fil}(angr).Parameter(lblind).data);
    dimm = parameterGr(angr).Parameter(lblind).dim;
    dimm(2) = length(parameterGr(angr).Parameter(lblind).data);
    parameterGr(angr).Parameter(lblind).dim = dimm;

    parameterGr(angr).Parameter(dind).data = ...
	cat(2,parameterGr(angr).Parameter(dind).data,...
	    parameters{fil}(angr).Parameter(dind).data);
    dimm = parameterGr(angr).Parameter(dind).dim;
    dimm(2) = length(parameterGr(angr).Parameter(dind).data);
    parameterGr(angr).Parameter(dind).dim = dimm;

    parameterGr(angr).Parameter(scind).data = ...
	cat(1,parameterGr(angr).Parameter(scind).data,...
	    parameters{fil}(angr).Parameter(scind).data);
    parameterGr(angr).Parameter(scind).dim = chind-1;

    parameterGr(angr).Parameter(ofind).data = ...
	cat(1,parameterGr(angr).Parameter(ofind).data,...
	    parameters{fil}(angr).Parameter(ofind).data);
    parameterGr(angr).Parameter(ofind).dim = chind-1;

    if fpgr
      parameterGr(fpgr).Parameter(cornind).data = ...
	  cat(3,parameterGr(fpgr).Parameter(cornind).data,...
	      parameters{fil}(fpgr).Parameter(cornind).data);
      parameterGr(fpgr).Parameter(cornind).dim = ...
	  size(parameterGr(fpgr).Parameter(cornind).data);

      parameterGr(fpgr).Parameter(orind).data = ...
	  cat(2,parameterGr(fpgr).Parameter(orind).data,...
	      parameters{fil}(fpgr).Parameter(orind).data);
      parameterGr(fpgr).Parameter(orind).dim = ...
	  size(parameterGr(fpgr).Parameter(orind).data);

      parameterGr(fpgr).Parameter(fpchind).data = ...
	  cat(2,parameterGr(fpgr).Parameter(fpchind).data,...
	      chind - 1-size(analogs{fil},2)...
	      +parameters{fil}(fpgr).Parameter(fpchind).data);
      parameterGr(fpgr).Parameter(fpchind).dim = ...
	  size(parameterGr(fpgr).Parameter(fpchind).data);

      parameterGr(fpgr).Parameter(zrind).data = ...
	  cat(2,parameterGr(fpgr).Parameter(zrind).data,...
	      parameters{fil}(fpgr).Parameter(zrind).data);
      parameterGr(fpgr).Parameter(zrind).dim = ...
	  size(parameterGr(fpgr).Parameter(zrind).data);

      parameterGr(fpgr).Parameter(tpind).data = ...
	  cat(1,parameterGr(fpgr).Parameter(tpind).data,...
	      parameters{fil}(fpgr).Parameter(tpind).data);
      parameterGr(fpgr).Parameter(tpind).dim = ...
	  size(parameterGr(fpgr).Parameter(tpind).data);

    end
  end
end

mrks = zeros(size(analog,1)*analogFrRt/vidFrRt,0,3);

% Write the file
writeC3D(mrks, vidFrRt, analog, analogFrRt, events, parameterGr,...
	 caminfo,residerr,file_out);


