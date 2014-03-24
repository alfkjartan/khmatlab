function vi = integrate_trapezoid(v, fs)
%
%  vi = integrate_trapezoid(v, fs)
%
% Integrates v using the trapezoid method. NaNs are taken to mean
% missing values, for which the width of the trapezoid is enlarged.
%
% Input
%    v   ->  vector
%    fs  ->  sampling frequency
% Output
%    vi  <-  integrated signal

% Kjartan Halvorsen
% 2008-04-15

vi = zeros(size(v));

for col = 1:size(v,2)
  for rad=2:size(v,1)
    vi(rad,col) = vi(rad-1, col) + ...
	
    
  end
end
