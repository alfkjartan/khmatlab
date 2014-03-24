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
%    numdef     ->  string on usual sprintf format for
%                     representing a num
%    num        ->  scalar.

% Kjartan Halvorsen
% 2003-10-17
  
  
if ( (column > length(tabledef)) | (column < 1) )
  error('Column index exceeds table limit');
end


colwidth = tabledef{column}{1};
colalign = tabledef{column}{2};
colsep = tabledef{column}{3};

s = sprintf(numdef, num);
if (length(s) > colwidth)
  warning('Number representation too wide for column');
  s(colwidth+1:end) = [];
end

margins = colwidth-length(s);

switch lower(colalign)
 case 'r'
  s = [blanks(margins), s];
 case 'c'
  lm = ceil(margins/2);
  rm = floor(margins/2);
  s = [blanks(lm), s, blanks(rm)];
 otherwise  % left align
  s = [s, blanks(margins)];
end

s = [s, colsep];