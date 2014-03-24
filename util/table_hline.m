function s = table_hline(tabledf)
  % s = table_hline(tabledef)
  % Returns a string with a horisontal line, the length being the
  % the same as the width of the whole table.
    
  tablewidth = 0;
  for c=1:length(tabledf)
    tablewidth = tablewidth + tabledf{c}{1};
    tablewidth = tablewidth + length(tabledf{c}{3});
  end
  
  s=chars(tablewidth, '-');
  