function s = column_string(tabledef, column, strng, colalign, columnspan)
% s = column_string(tabledef, column, string, align, columnspan)
% Returns a string with the given double in the format for the
% given column. Padded with spaces as needed.
%
% Type help columndef to see information on the table definition
%
% Input
%    tabledef     ->  Cell array describing the table format
%    column       ->  the column index
%    string       ->  string to print
%    align        ->  optional. Will override definition in tabledef.
%    columnspan   ->  optional. number of columns that the string
%                     should span

% Kjartan Halvorsen
% 2003-10-17
  
  
if ( (column > length(tabledef)) | (column < 1) )
  error('Column index exceeds table limit');
end


colwidth = tabledef{column}{1};
if ((nargin < 4) | isempty(colalign))
  colalign = tabledef{column}{2};
end

colsep = tabledef{column}{3};

s=strng;

if ( (nargin == 5) & (columnspan > 1) )
  if ( column + columnspan -1 <= length(tabledef))
    for c=(column+1):(column+columnspan-1)
      colwidth = colwidth + tabledef{c}{1} + length(tabledef{c}{3});
    end
    colsep = tabledef{c}{3};
  end
end

if (length(s) > colwidth)
  warning('String too wide for column');
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