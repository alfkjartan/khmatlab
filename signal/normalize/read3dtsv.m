function [attr,data]=read3dtsv(fid)
% [attr,data]=read3dtsv(fid)
% Function that reads a tsv file with 3D marker data. 



% Read the first attribut name and value
line=fgetl(fid);

attr={};
acnt=0;

[attrname,attrval]=strtok(line);

while ((line ~= -1) & isempty(str2num(attrname))) 
   % first line does not start with number and not EOF
   if ( ~isempty(attrval) )
     attrval(1)=[];
     attrval=deblank(attrval);
     attrval=fliplr(deblank(fliplr(attrval)));
   else
     attrval = 'NA';
   end
   
   acnt=acnt+1;
   attr=putvalue(attr,attrname,attrval);
   line=fgetl(fid);
   [attrname,attrval]=strtok(line);  
end

if (feof(fid))
   error('Unexpected EOF');
end

% Parse the marker names to construct a cell array with a string for
% each marker
nmarkers=sscanf(getvalue(attr,'NO_OF_MARKERS'),'%d');
mnstr=getvalue(attr,'MARKER_NAMES');

if (~isempty(mnstr))
  mnames=cell(nmarkers,1);
  for m=1:nmarkers
    [mnames{m},mnstr]=strtok(mnstr);
  end
  attr=putvalue(attr,'MARKER_NAMES',mnames); 
end


% Parse the channel names (if exists) to construct a cell array
% with a string for each channel

nchannels=sscanf(getvalue(attr,'TOT_NO_OF_CHAN'),'%d');
chnstr=getvalue(attr,'CHANNEL_NAMES');

if (~isempty(chnstr))
  chnames=cell(nchannels,1);
  for m=1:nchannels
    [chnames{m},chnstr]=strtok(chnstr);
  end
  attr=putvalue(attr,'CHANNEL_NAMES',chnames); 
end

%disp('Read attributes')

nframes=sscanf(getvalue(attr,'NO_OF_FRAMES'),'%d');

data=sscanf(fixline(line),'%f');
l=1;
while 1
  line = fgetl(fid);
  if ~ischar(line), break, end
  fline = fixline(line);
  if isempty(fline), break, end
  try
    data = cat(2, data, sscanf(fline,'%f'));
  catch
    error(['Unexpected number of elemnts on line: ',int2str(l)]);
  end
  l=l+1;
end

data = data';

function nl = fixline(l)
  % Replaces commas with periods.
  pat = ',';
  nl = regexprep(l,pat,'.');
  
