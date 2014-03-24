function v = centraldiff(x, fs)
%  v = centraldiff(x, fs)
%  Computes velocity estimates as the central difference of x.
%  Operates along the columns of x. v has the same size as x. 
%  The first and last elements are estimated by the one-step ahead
%  and one-step back approximation. fs is the sampling frequency.
%

% Kjartan Halvorsen
% 2004-03-23

v=nan(size(x));
for d=1:size(x,3)
  xplus1 = x(3:end,:, d);
  xmin1 = x(1:end-2,:, d);

  v(1,:,d) = (x(2,:,d) - x(1,:,d))*fs;
  v(2:(end-1),:, d) =  (xplus1 - xmin1)*fs*0.5;
  v(end, :, d) = (x(end,:,d) - x(end-1,:,d))*fs;
end
