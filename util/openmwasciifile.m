function [data, dt, filename] = openmwasciifile(varargin)
% Opens a file dialog, reads the content of the text file, assuming
% it is exported from power lab.
  
% Kjartan Halvorsen
% 2003-12-02
  
if (nargin == 0)
  title = 'Pick a data file';
else
  title = varargin{1};
end

  [filename, pathname] = uigetfile('*.*', title);

  if (~filename) 
    data=[];
    dt = [];
    return; 
  end

  fid = fopen(fullfile(pathname, filename));

  dt=[];
  tab = sprintf('\t');
  while 1
    headerline = fgetl(fid);
    try 
      [tok, rest] = strtok(headerline, ['= ', tab]);
      if (strfind(tok, 'SamplingFreq')) % Sampling interval
	dt = str2num(strtok(rest));
      end
      if ~isnan(str2double(tok))
	break
      end
    end
  end
  
  firstline = sscanf(headerline, '%f\t%f\t%f');
  ncols = length(firstline);
  data = fscanf(fid, '%f');
  fclose(fid);

  data = reshape(data, ncols, length(data)/ncols);
  data = cat(1, firstline', data');

