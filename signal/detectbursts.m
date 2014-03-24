function [brstfr]=detectbursts(md, threshold)
%  [brstfr]=detectbursts(md, threshold)
% 
% Detect bursts in emg data. Returns the sample indices where the
% local variance in the signal exceeds <threshold * baseline
% variance>.
%
% Input
%    md              ->   emg data (N x 1), recitified and low pass
%                         filtered
%    threshold       ->   Threshold for detection of burst.

% Kjartan Halvorsen
% 2003-09-11
%
% Revisions
% 2004-10-15  Assume emg data already rectified and low pass filtered

% Parameter that decides how many samples, in % of the cycle
% frequency, that must be above the threshold for the burst to be detected
burstwindow = 5;

nfr = length(md);

% Find periodicity in the signal.

% Use the peridogram to estimate the cycle frequency.
mdfft=fft(md,48*nfr);
psd=abs(mdfft(1:24*nfr));

maxmdf=max(psd(48:end));
peak=find(psd==maxmdf);

%plot(psd)

%keyboard

T=ceil(48*nfr/(peak(1)-1));

% The dominating frequency is found above. Use this to generate a
% cosine wave of the same length as the signal. Multiply
% with the signal and take the mean. The result, vv will be a
% function of the phase of the signal.  If we define the
% phase to be zero for a sine formed signal, then the
% phase will be ph = asin(vv*2);

sw=cos(linspace(0,2*pi*(peak(1)-1)/48,length(md)));

prod = max( min( 2*mean(sw.*md')/max(md), 1), -1);

ph=fix(asin(prod)/(2*pi)*T);

peaks=(1:T:(1.2*nfr))'+ph;

while ((peaks(end)>nfr) | peaks(1)<1)
  if (peaks(end)>nfr) 
    if peaks(end)>nfr+T/2
      peaks(end)=[];
    else
      peaks(end)=nfr;
    end
  elseif (peaks(1)<1)
    if peaks(1)<-T/2
      peaks(1)=[];
    else
      peaks(1)=1;
    end
  end
end

% Find the start and end of each burst.

% Find first the rms of the baseline part
baseline = clipmdata(md(1:min(fix(8*T),end)),...
		     'Indicate region with baseline signal');
rms_base = mean(baseline);
thr = rms_base * threshold;

%brstfr = cell(0);
brstfr = [];

bfrs = floor(burstwindow*T/100);

for i=1:length(peaks)-1
  startind=0;
  endind=0;
  k=-1;
  while((startind == 0) & (peaks(i) + k < peaks(i+1)))
    k = k+1;
    if (md( (max(1, peaks(i)+ k - bfrs)) : peaks(i)+ k) > thr)
      startind = max(1, peaks(i) + k - bfrs);
    end
  end
  while((endind == 0) & (peaks(i) + k < peaks(i+1)))
    k = k+1;
    if (md( (max(1, peaks(i)+ k - bfrs)) : peaks(i)+ k) < thr)
      endind = max(1, peaks(i) + k - bfrs);
    end
  end
  %brstfr = cat(1, brstfr,  {[startind endind]});
  brstfr = cat(1, brstfr,  startind, endind);
end

