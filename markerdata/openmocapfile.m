function [mdata,adata,fname] = openmocapfile(varargin)
% Convenient way of loading a file with mocap data. The file format
% can either be tsv or c3d.
%
% Usage:
%     [mdata, adata, fname] = openmocapfile
%     [mdata, adata, fname] = openmocapfile(titlestring)
%     [mdata, adata] = openmocapfile(titlestring, fname)
%
% Based on opentsv, and returns data in the same format.
  
% Kjartan Halvorsen
% 2004-07-29

  if (nargin == 0)
    prompt = 'Pick a motion capture data file';
  else
    prompt = varargin{1};
  end
  
  if (nargin == 2)
    fname = varargin{2};
  else
    [filename, pathname] = uigetfile(...
	{'*.tsv;*.TSV;*.c3d;*C3D', 'MoCap files (*.tsv, *.c3d)'}, prompt);
    if (~filename) 
      mdata=[];
      adata=[];
      fname=[];
      return; 
    end

    fname = fullfile(pathname, filename);
  end

  [pthpart, namepart, ext] = fileparts(fname);

  if ( strcmp(ext,'.tsv') | strcmp(ext,'.TSV') | strcmp(ext, '.txt'))
    % Same as before in opentsv.m
    fid = fopen(fname);
  
    [attr, md] = read3dtsv(fid);
    mdata = {attr, md/1000}; % Divide by 1000 to convert from mm to m
    fclose(fid);

    if (nargout > 1)
      % Try to open the file with analog data.
      fid = fopen( fullfile(pthpart, [namepart, '_a', ext]) );
      
      if (fid > 0)
	[attra, ad] = read3dtsv(fid);
	adata = {attra, ad};
	fclose(fid);
      else
	adata = {};
      end
    end
  
  elseif ( strcmp(ext,'.c3d') | strcmp(ext,'.C3D') )
    try
      [markers, vidfrate, analogs, analogfrate, evnt, ParameterGroup] = ...
	  readC3D(fname);
    catch
      mdata=[];
      adata=[];
      warning(['File not found: ', fname])
      return
    end
    
    % Convert to the same format as returned by earlier opentsv.
    
    % Marker data
    
    % The attributes. Only the ones I know are used
    [slsk,pointgrind] = intersect([ParameterGroup.name], 'POINT');
    [slsk,lblsind] = intersect([ParameterGroup(pointgrind).Parameter.name],...
			       'LABELS');
    markattr = putvalue('NO_OF_FRAMES',num2str(size(markers,1)));
    markattr = putvalue(markattr,'NO_OF_MARKERS',num2str(size(markers,2)));
    markattr = putvalue(markattr,'FREQUENCY',num2str(vidfrate));

    if isempty(lblsind)
      markattr = putvalue(markattr,'MARKER_NAMES',{ });
    else
      markattr = putvalue(markattr,'MARKER_NAMES',...
			  ParameterGroup(pointgrind).Parameter(lblsind).data');
    end
      
    mdata = {markattr, ...
	     reshape(permute(markers,[1 3 2]), ...
		     [size(markers,1) size(markers,2)*3])};
    
    
    % Analog data
    
    % The at1tributes. Only the ones I know are used
    anattr = putvalue('NO_OF_SAMPLES',num2str(size(analogs,1)));
    anattr = putvalue(anattr,'TOT_NO_OF_CHAN',num2str(size(analogs,2)));
    anattr = putvalue(anattr,'FREQUENCY',num2str(analogfrate));
    anattr = putvalue(anattr,'CHANNEL_NAMES',...
			ParameterGroup(2).Parameter(3).data');

    adata = {anattr, analogs};
  
  else
    mdata=[];
    adata=[];
    warning(['Unknown extension: ', ext])
  end
  
  
  % Check marker names, add numbers if nonunique
  
  if ~isempty(mdata)
    mnames = getvalue(mdata{1}, 'MARKER_NAMES');
    uniquenames = unique(mnames);
    
    if ( length(uniquenames) < length(mnames) )
      for i=1:length(uniquenames)
	ind = find(strcmp(mnames, uniquenames{i}));
	if length(ind)>1
	  for j=1:length(ind)
	    mnames{ind(j)} = strcat(mnames{ind(j)},int2str(j));
	  end
	end
      end
      mdata{1} = putvalue(mdata{1}, 'MARKER_NAMES', mnames);
    end

  end
  
  
