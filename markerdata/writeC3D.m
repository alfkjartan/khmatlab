function byteswritten = writeC3D(Markers,VideoFrameRate, ...
				 AnalogSignals, AnalogFrameRate, ...
				 Event,ParameterGroup,CameraInfo, ...
				 ResidualError, FullFileName)
%  byteswritten = writeC3D(Markers,VideoFrameRate, ...
%				 AnalogSignals, AnalogFrameRate, ...
%				 Event,ParameterGroup,CameraInfo, ...
%				 ResidualError, FullFileName)
%
% Writing 3d data to a c3d file. Based on readC3D by Alan Morris
% and Jaap Haarlar
%
% Input
%    Markers         ->  3D marker data (Nvideoframes x Nmarkers x 3)
%    VideoFrameRate  ->  Frames/sec
%    AnalogSignals   ->  Analog data (Nanalogsamples x Nsignals)
%    AnalogFrameRate ->  Samples/sec
%    Event           ->  Struct. Event(Nevents).time ..value ..name
%    ParameterGroup  ->  Struct. 
%                        ParameterGroup(Ngroups).Parameters(Nparameters)
%                                                         .data ..etc
%    CameraInfo      ->  Marker related info (NvideoFrames x Nmarkers)
%    ResidualError   ->  Marker related residual info (NvideoFrames x Nmarkers)
%    FullFileName    ->  file (including path) to write to 
%

% Kjartan Halvorsen
% 2004-08-17

% The processor type 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)
proctype = 1;

%-------------------------------------------------------------
% 
%     open the file
% 
%-------------------------------------------------------------

  ind=findstr(FullFileName,'\');
  if ind>0 
    FileName=FullFileName(ind(length(ind))+1:length(FullFileName)); 
  else 
    FileName=FullFileName; 
  end

  fid=fopen(FullFileName,'wb','n'); % native format (PC-intel)

  if fid==-1,
    h=errordlg(['File: ',FullFileName,' could not be opened'],'application error');
    uiwait(h)
    return
  end

%-------------------------------------------------------------
% 
%     write header
% 
%-------------------------------------------------------------

  % Write the record number of the  parameter section. Always 2
  NrecordFirstParameterblock=2;
  fwrite(fid, NrecordFirstParameterblock, 'int8');     
  key = 80;
  fwrite(fid, key, 'int8');                           % key = 80;

  % Number of markers
  fwrite(fid, size(Markers,2), 'int16');
  
  % Number of analog measurements per frame
  fwrite(fid, AnalogFrameRate / VideoFrameRate * size(AnalogSignals,2), ...
	 'int16');
  
  % First frame (always 1) and last of marker data.
  firstField = 1;
  lastField = size(Markers, 1);
  fwrite(fid, firstField, 'int16');
  fwrite(fid, lastField, 'int16');
  
  % Max interpolation gap. Default 10.
  maxgap = 10;
  fwrite(fid, maxgap, 'int16');
  
  % Scaling factor
  Scale = -0.08;
  fwrite(fid, Scale, 'float32');
  
  % Before we can write the starting record number of the marker
  % and analog data, we need to know how many records the parameter
  % section needs.
  tmpParamfile = ['tmpParameterFile',date,'.c3d'];
  fid2 = fopen(tmpParamfile,'wb');
  [byteswritten, parameterblocks] = writeParams(fid2,ParameterGroup)
  fclose(fid2);

  firstdatarecord = 2 + parameterblocks;
  
  %fseek(fid, 16, 'bof');
  fwrite(fid, firstdatarecord, 'int16');
  
  % Number of analog samples per frame
  NanalogFramesPerVideoFrame = AnalogFrameRate / VideoFrameRate;
  fwrite(fid, NanalogFramesPerVideoFrame, 'int16');
  
  % Sampling frequency of marker data
  fwrite(fid, VideoFrameRate, 'float32');
  
  % write zeros, empty part of header
  fwrite(fid, zeros(298-ftell(fid),1), 'int8');
  
  
  %-------------------------------------------------------------
  % 
  %     Write events
  % 
  %-------------------------------------------------------------
  
  if isempty(Event)
    fwrite(fid, 0, 'int16'); % Event indicator
  else
    fwrite(fid, 12345, 'int16');
    fwrite(fid, length(Event), 'int16');
    fwrite(fid, 0, 'int16'); % Reserved for future use
    for i = 1:length(Event)
      fwrite(fid, Event(i).time, 'float32');
    end
    
    % Write junk up until event value part
    fwrite(fid, zeros(188*2 - ftell(fid),1), 'int8');
    for i = 1:length(Event)
      fwrite(fid, Event(i).value, 'int8');
    end

    % Write junk up until event name part
    fwrite(fid, zeros(198*2 - ftell(fid),1), 'int8');
    for i = 1:length(Event)
      fwrite(fid, double(Event(i).name{1}), 'char');
    end
  end

  fwrite(fid, zeros(512-ftell(fid), 1), 'int8');
  
  %-------------------------------------------------------------
  % 
  %     Write / copy parameters
  % 
  %-------------------------------------------------------------
  
  %keyboard
  % Copy the parameters from the temporary parameter file 

  % Write the header part of the parameter section
  if (ftell(fid) ~= 512)
    warning(['Expected file pointer to be at 512, but it is ', ...
	     int2str(ftell(fid))])
  end
  
  fwrite(fid, 0, 'int8');
  fwrite(fid, 80, 'int8'); % Key. Must be 80
  fwrite(fid, parameterblocks, 'int8'); 
  fwrite(fid, 84, 'int8'); % 83+ processor type. 1: intel

  % Copy the rest
  fid2 = fopen(tmpParamfile,'rb');

  pblock = fread(fid2, 512*parameterblocks-4, 'int8');
  fwrite(fid, pblock, 'int8');
    
  fclose(fid2);
  
  

  %-------------------------------------------------------------
  % 
  %     Write data
  % 
  %-------------------------------------------------------------


  posatstartdata = ftell(fid);
  firstdatarecord = firstdatarecord;
  firstdatawritepos = (firstdatarecord-1)*512;
  sizepblock = size(pblock);
  parameterblocks = parameterblocks;
  
  NvideoFrames = size(Markers, 1);
  Nmarkers = size(Markers, 2);
  NanalogChannels = size(AnalogSignals, 2);

  dataPerFrame =   4*Nmarkers + NanalogChannels*NanalogFramesPerVideoFrame ;

  nodataind = find(ResidualError==-1);
  reserr = ResidualError/abs(Scale);
  reserr(nodataind) = 255;
  CameraInfo(nodataind) = 255;
  CamRes = ...
      double(reshape( ...
	  (CameraInfo*256 + reserr)',...
	  1, Nmarkers*NvideoFrames));
  Markers = cat(1,...
		reshape(permute(Markers, [3 2 1]),...
			3, Nmarkers*NvideoFrames) ,...
		CamRes);
  AllData = cat(1, reshape( Markers, 4*Nmarkers, NvideoFrames),...
		reshape(AnalogSignals',...
			NanalogChannels*NanalogFramesPerVideoFrame,...
			NvideoFrames));
  
  fwrite(fid,AllData,'float32');
  
  fclose(fid)
    
  
  
  
  
function [byteswritten, blockswritten] = writeParams(fid, ParameterGroup)
  % Write parameter section header
% fseek(fid, 512, 'bof');
%  fwrite(fid, 0, 'int8');
%  fwrite(fid, 80, 'int8'); % Key. Must be 80
%  fwrite(fid, 0, 'int8'); % Will later hold the number of parameter
                          % blocks
%  fwrite(fid, 84, 'int8'); % 83+ processor type. 1: intel
  
  % write the groups and parameters
  byteswritten = 4;
  %byteswritten = 0;
  
  for i = 1:length(ParameterGroup)
    byteswritten = byteswritten + writeGroup(fid, ...
					     ParameterGroup(i), ...
					     -i);
  end
  
  % Blocks used to write the parameter section
  blockswritten = ceil(byteswritten / 512);
  ftell(fid)
  

  % fseek(fid, 2, 'bof');
 % fwrite(fid, blockswritten, 'int8');

  % Clear the remainder of the last 512 byte block
  %fseek(fid, byteswritten, 'bof');
  fwrite(fid, zeros((blockswritten)*512-byteswritten, 1), 'int8');
  
  ftell(fid)
  
  
  
function byteswritten = writeGroup(fid, Group, id)
  namelength = length(Group.name{1}); 
  Group.description{1} = fixchars(Group.description{1});
  descrlength = length(Group.description{1}); 
  fwrite(fid, namelength, 'int8');
  fwrite(fid, id, 'int8');
  fwrite(fid, double(Group.name{1}), 'char'); 
  fwrite(fid, descrlength+3, 'int16') % The offset
  fwrite(fid, descrlength, 'int8');
  fwrite(fid, double(Group.description{1}), 'char'); 

  byteswritten = 5 + namelength + descrlength;
  
  for i = 1:length(Group.Parameter)
    byteswritten = byteswritten + writeParameter(fid, Group.Parameter(i), ...
						 abs(id));
  end
  
      
  
function byteswritten = writeParameter(fid, Parameter, groupID)

  if (isfield(Parameter, 'description'))
    if iscell(Parameter.description)
      Parameter.description{1} = fixchars(Parameter.description{1});
    else
      Parameter.description{1} = fixchars(Parameter.description);
    end
  else
    Paremeter.description = {'no description'};
  end
  
  if ~isfield(Parameter,'dim')
    Parameter.dim=1;
  end
  if isempty(Parameter.dim)
    Parameter.dim = 1;
  end
  
  namelength = length(Parameter.name{1}); 
  descrlength = length(Parameter.description{1}); 
  
  % Determine length of data
  dlen = abs(Parameter.datatype);
  for d=1:length(Parameter.dim)
    dlen = dlen*Parameter.dim(d);
  end
  
  fwrite(fid, namelength, 'int8');
  fwrite(fid, groupID, 'int8');
  fwrite(fid, double(Parameter.name{1}), 'char'); 

  % The offset:
  offset = 5 + length(Parameter.dim) + dlen + descrlength;
  fwrite(fid, offset, 'int16');
  %keyboard
  
  fwrite(fid, Parameter.datatype, 'int8');
  fwrite(fid, length(Parameter.dim), 'int8');
  fwrite(fid, Parameter.dim, 'int8');

  %if (strcmp(Parameter.name,'ZERO'))
  %  keyboard
  %end
  
  % Write data
  switch Parameter.datatype
   case -1  % Char
    % TODO 2004-10-12:
    % Fix so that truncated character arrays do fill up the number
    % of characters that Parameter.dim specify.
    % Probably safe to assume that possible dimensions are 1 or 2.
    wordlength = Parameter.dim(1);
    if (length(Parameter.dim) == 2)
      for i = 1:length(Parameter.data)
	fwrite(fid, cat(2, double(Parameter.data{i}), ...
			repmat(32,1, ...
			       wordlength-length(Parameter.data{i}))), 'char');
      end
    else
      Parameter.name
      Parameter.data
      fwrite(fid, cat(2, double(Parameter.data{1}), ...
		      repmat(32,1, ...
			     wordlength-length(Parameter.data{1}))), 'char');
    end
   case 1  % boolean
    fwrite(fid, Parameter.data, 'int8');
   case 2  % integer
    fwrite(fid, Parameter.data, 'int16');
   case 4 % Float
    fwrite(fid, Parameter.data, 'float32');
  end
  
  fwrite(fid, descrlength, 'int8');
  fwrite(fid, double(Parameter.description{1}), 'char'); 

  byteswritten = offset + 2 + namelength;
  
  %keyboard  
  
function strout = fixchars(strin)
  % Removes characters that are not ascii
  strdbl = double(strin);
  strdbl(find(strdbl>126)) = [];
  strout = char(strdbl);
  
