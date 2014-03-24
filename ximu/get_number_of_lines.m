function nl = get_number_of_lines(filestr)
  %% Returns the number of lines in the file

  [status, wcstr] = system(["wc -l ", filestr]);
  [nl] = sscanf (wcstr, "%f ", 1);
