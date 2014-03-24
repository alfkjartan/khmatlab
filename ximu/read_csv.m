function [dta, header] = read_csv(filestr)
  %% Reads a comma separated text file, assuming the first line is a
  %% head line and data starting at line 2.
  %%
  %% Input
  %%    filestr      ->  string with file name
  %% Output
  %%    dta          ->  data in matrix
  %%    header       ->  cell array with column description

  %% Kjartan Halvorsen
  %% 2012-09-11


  %% Get number of lines to read, to speed things up
  nlines = get_number_of_lines(filestr);

  %% Read the first line. Parse the header
  fid = fopen(filestr, 'r');
  line1 = fgetl(fid);
  header = strsplit(line1, ",");
  
  ncolumns = length(header);

  template = strcat(repmat("%f,", 1, ncolumns-1), "%f\n");

  %% The total number of data elements in the file
  nelements = (nlines-1) * ncolumns;

  %% Read the data
  [vals, count] = fscanf(fid, template, nelements);

  if (count != nelements)
    warning("Unexpected number of values read!")
  end

  dta = (reshape(vals, ncolumns, count/ncolumns))';

