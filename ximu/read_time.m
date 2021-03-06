function [packet_number_time] = read_time(filestr)
  %% Reads the file, assuming it is a DateTime file exported by
  %% ximu-gui. Returns the packet number and corresponding time as the
  %% number of seconds since the epoch.

  %% Kjartan Halvorsen
  %% 2012-09-10

  %% The file contains 7 columns:
  %% packet,  year,  month,  day,  hour,  minute,  second

  %% First, get number of lines to read, to speed things up

  nlines = get_number_of_lines(filestr);

  dttme = dlmread(filestr, ',', [1, 0, nlines-1, 6]);

  packet_number_time = zeros(size(dttme, 1), 2);
  
  packet_number_time(:,1) = dttme(:,1);
  
  for i=1:size(dttme,1)
    packet_number_time(i,1) = datenum(dttme(i,2), \
				      dttme(i,3), \
				      dttme(i,4), \
				      dttme(i,5), \
				      dttme(i,6), \
				      dttme(i,7));
  end

