function s = column_num(tabledef, column, numdef, num)
% s = column_num(columndef, column, numdef, num)
% Returns a string with the given double in the format for the
% given column. Padded with spaces as needed.
%
% Type help columndef to see information on the table definition
%
% Input
%    tabledef     ->  Cell array describing the table format
%    column       ->  the column index
%    floatdef     ->  string on usual sprintf format for
%                     representing a float
%    float        ->  scalar.

% Kjartan Halvorsen
% 2003-10-17
  
  
if ( (column > length(tabledef)) | (column < 1) )
  error('Column index exceeds table limit');
end


colwidth = tabledef{column}{1};
colalign = tabledef{column}{2};
colsep = tabledef{column}{3};

s = sprintf(floatdef, float);
if (length(s) > colwidth)
  warning('Float representation too wide for column');
  s(colwidth+1:end) = [];
end

margins = colwidth-length(s);

switch lower(colalign)
 case 'r'
  s = [spaces(margins), s];
 case 'c'
  lm = ceil(margins/2);
  rm = floor(margins/2);
  s = [spaces(lm), s, spaces(rm)];
 otherwise  % left align
  s = [s, spaces(margins)];
end

s = [s, colsep];