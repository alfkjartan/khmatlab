function export_values(val, allnames, names2export, fname, frame2export)
%  export_values(val, allnames, names2export, fname, frame2export)
% Exports the rows and the column of val indicated by the names and
% frame number.
% The exported text file will have as many columns as there are
% names2export (pluss one for the label of the event) and as many
% rows as are indicated by frame2export. 
%
% Input
%    val        ->  The values (nst x nfrs)
%    allnames   ->  The names, a (nst x 1) cell array of strings
%    names2export  ->  Cell array with names of variables to export.
%    fname      ->  The name of the text file to write to
%    frame2export  ->  Two possibilities: vector of indices or cell
%                      array where the first column contains names
%                      of the event, and the second the
%                      corresponding index (column index to val).

% Kjartan Halvorsen
% 2009-07-08
