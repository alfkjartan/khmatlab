function [mdata, adata, fname] = opentsv(varargin)
% Convenient way of loading a file with tsv marker data
% Usage:
%     [mdata, adata, fname] = opentsv
%     [mdata, adata, fname] = opentsv(titlestring)
%     [mdata, adata] = opentsv(titlestring, fname)
%
% OBSOLETE AS of 2004-07-29. USE openmocapfile INSTEAD.

% Kjartan Halvorsen
% 2003-11-06
    
  if (nargin == 0)
    prompt = 'Pick a tsv data file';
  else
    prompt = varargin{1};
  end
  
  if (nargin == 2)
    fname = varargin{2};
  else
    [filename, pathname] = uigetfile('*.*', prompt);
    if (~filename) 
      mdata=[];
      adata=[];
      fname=[];
      return; 
    end

    fname = fullfile(pathname, filename);
  end


  fid = fopen(fname);
  
  [attr, md] = read3dtsv(fid);
  mdata = {attr, md};
  fclose(fid);

  if (nargout > 1)
    % Try to open the file with analog data.
    [pthpart, namepart, ext] = fileparts(fname);
    fid = fopen( fullfile(pthpart, [namepart, '_a', ext]) );

    if (fid > 0)
      [attra, ad] = read3dtsv(fid);
      adata = {attra, ad};
      fclose(fid);
    else
      adata = {};
    end
  end
  
  