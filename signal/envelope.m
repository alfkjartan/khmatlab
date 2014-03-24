function de=envelope(d, window, type)
%  de=envelope(d, param, [type])
% Computes the envelope of the signal as a lowpass filtered version
% of the rectified signal, or, if type == 'rms' as the local root-mean-square
% of the signal in a window of length given, centered on each data
% point.

% Kjartan Halvorsen
% 2003-09-11
% Revisions
% 2004-10-15  Changed to lowpass filtering of rectified signal
%             if third argument is given, must be 'rms' to get old
%             results with rms windowing.

filterorder = 2;

if (nargin < 3)
  type = 'lp';
end


[N,cols] = size(d);

if strcmp(type, 'rms')
  de = zeros(size(d));
  wl = fix(window/2);

  divisor = sqrt(2*wl+1);
  for c=1:cols
    for i=1:N
      wind = d(max(1,i-wl):min(N,i+wl),c);
      de(i,c) = norm(wind)/sqrt(length(wind));
    end
  end
else
  [b,a] = butter(filterorder, window);
  de = filtfilt(b,a, abs(d));
end

  