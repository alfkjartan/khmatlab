function [mmd, ind, nfound] = extractmarkers(varargin)
% [mmd, ind, nfound]=extractmarkers(md,mnames,submnames) or
% [...]=extractmarkers(mdata, submnames)
% Returns a matric containing the set of columns of md that corresponds
% to the cell array of marker names submnames.
% If the cell array mnames is not given, then the first argument
% must be a cell array with attributes and marker data: mdata =
% {attr, md}.
%
% If marker names in submnames are missing from mnames, the
% corresponding columns in the returned matrix contain zeros
% Output
%   mmd    ->   matrix with marker data
%   ind    ->   indices into original marker data
%   nfound ->   the number of markers found.

%
% Kjartan Halvorsen
% 2000-11-08
%
% Revisions
% 2003-08-25   Changed to deal with repeated markernames in submnames 
% 2003-11-17   Can take two arguments, first being {attr, md}.
% 2003-12-02   Can take multiple mdata arrays
% 2005-04-25   Change to be caseinsensitive.
% 2008-01-09   Changed to return zeros when marker not found in data
% 2008-07-07   Changed to return also the number of markers found

if ( iscell(varargin{1}) & (length(varargin{1}) > 1) )
  if ( iscell(varargin{1}{2}) ) % Set of mdata 
    mmd=[];
    for i=1:length(varargin{1})
      md = extractmarkers(varargin{1}{i}, varargin{2:end});
      mmd = cat(1, mmd, md);
    end
    return
  end
end

if ( nargin == 2 & iscell(varargin{1}) )
  md = varargin{1}{2};
  mnames = getvalue(varargin{1}{1}, 'MARKER_NAMES');
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

mnames = upper(mnames);
submnames = upper(submnames);

% Make a copy with three extra columns containing zeros. When a
% requested marker is missing from the data set, the corresponding
% index will be set to the last three columns, thus returning zeros
% for that marker. This guarantees that the matrix returned has
% a predictable size.

ncols = size(md, 2);
mdcopy = cat(2, md, zeros(size(md,1), 3));
ncols2 = size(mdcopy, 2);

ind = [];
nfound = 0;
for i = 1:length(submnames)
  [cmn,mind]=intersect(mnames,submnames(i));
  if isempty(cmn)
    ind=[ind ncols2-2:ncols2];
  else
    nfound = nfound+1;
    ind=[ind (3*(mind-1)+1):(3*mind)];
  end

end

mmd=mdcopy(:,ind);
