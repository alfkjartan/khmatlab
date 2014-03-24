function [adata] = openplianceasc(fname)
% Loads a text file with data from the pliance system.
  
% Kjartan Halvorsen
% 2012-02-22

fid = fopen(fname);
  
%% Read header. Look for "total time [secs]" and "time per frame [secs]"
%% to get the number of lines to read in.

line=fgetl(fid);

[frst,rem]=strtok(line);

tab = sprintf('\t');

total_time = 0;
time_per_frame = 0;
got_total_time = 0;
got_time_per_frame = 0;

headerlines = 0;
while ((line ~= -1) & isnan(str2double(frst))) 
  headerlines = headerlines + 1;
  cols = strsplit(line,tab);

  for c=cols
    if ~got_time_per_frame
      tok = regexp(c, '^time per frame \[secs\]:\s+(\d+\.\d+)$', 'tokens');
      if (~isempty(tok{1}))
	got_time_per_frame = 1;
	time_per_frame = str2double(tok{1}{1}{1});
      end
    end
    if ~got_total_time
      tok = regexp(c, '^total time \[secs\]:\s+(\d+\.\d+)$', 'tokens');
      if (~isempty(tok{1}))
	got_total_time = 1;
	total_time = str2double(tok{1}{1}{1});
      end
    end
  end

  line=fgetl(fid);
  [frst,rem]=strtok(line);
  while isempty(frst)
    headerlines = headerlines + 1;
    line = fgetl(fid);
    [frst,rem]=strtok(line);
  end
end
if (feof(fid))
   error('Unexpected EOF');
end

fclose(fid);

attr = {};
attr = putvalue(attr, 'FREQUENCY',  1.0/time_per_frame);


data = dlmread(fname, tab, headerlines, 0);

%% Reorder the columns to match our definition
ourorder = flipud(reshape(2:65, 4, 16));
ourorder = ourorder(:)';
adata = {attr, data(:,[1 ourorder])};
%%adata = {attr, data(:,1:65)};