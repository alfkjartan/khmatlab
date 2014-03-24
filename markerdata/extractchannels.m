function mmd=extractchannels(varargin)
% mmd=extractchannels(md,mnames,submnames) or
% mmd=extractchannels(mdata, submnames)
% Returns a matrix containing the set of columns of md that corresponds
% to the cell array of channel names submnames.
% If the cell array mnames is not given, then the first argument
% must be a cell array with attributes and marker data: mdata =
% {attr, md}.

% Based on extractmarkers
% Kjartan Halvorsen
% 2004-01-26
%

  
if ( iscell(varargin{1}) & (length(varargin{1}) > 1) & ...
     iscell(varargin{1}{2}) ) % Set of mdata 
  mmd=[];
  for i=1:length(varargin{1})
    md = extractchannels(varargin{1}{i}, varargin{2:end});
    mmd = cat(1, mmd, md);
  end
  return
end

if ( nargin == 2 & iscell(varargin{1}) )
  md = varargin{1}{2};
  mnames = getvalue(varargin{1}{1}, 'CHANNEL_NAMES');
  submnames = varargin{2};
else
  if iscell(varargin{1})
    md = varargin{1}{2};
  else
    md = varargin{1};
  end
  mnames = varargin{2};
  submnames = varargin{3};
end
  
if ischar(submnames)
  submnames = {submnames};
end

ind = [];
for i = 1:length(submnames)
  [cmn,mind]=intersect(mnames,submnames(i));
  ind=[ind mind];
end

mmd=md(:,ind);
