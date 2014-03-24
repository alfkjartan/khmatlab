function nmd=filtermdata(md,samplefreq,cutofffreq, order)
%  nmd=filtermdata(md,samplefreq,cutofffreq)
% Filters the data with a butterworth filter of given order at the given
% cutoff frequency.

% Kjartan Halvorsen
% 2000-09-19

nmd=detrend(md);
trend=md-nmd;

[b,a]=mybutter(order/2,2*cutofffreq/samplefreq);

nmd=filter(b,a,nmd);
nmd=flipud(nmd);
nmd=filter(b,a,nmd);
nmd=flipud(nmd);

nmd=nmd+trend;

