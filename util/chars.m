function s = chars(n, strng)
  % s = chars(n, ch)
  % Returns a string with a number of the given character,
  
  % Kjartan Halvorsen
  % 2003-10-17

  if ( nargin < 2)
    strng = ' ';
  end
  
  s='';
  for i=1:n
    s=[s, strng];
  end
  
  