function s = strtrim(str)
 % Trims white space from the beginning and end of the string.
   
 % Kjartan Halvorsen
 % 2003-12-02
   
 s = fliplr( deblank( fliplr( deblank(str) ) ) );
 