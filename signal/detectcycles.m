function [stpfr, mdvel, mdacc]=detectcycles(md, T, fs, useseries)
% [stpfr, dvel, dacc]=detectcycles(md, T, fs, useseries)
% 
% Detect cycles in the data md. 
% The data is first differentiated. The start
% of the cycle is taken to be at the lower peak of the velocity.
%
%
% Input
%    md              ->   displacement data (N x 1)
%    T               ->   Approximate period length (in # of samples)
%    fs              ->   sampling frequency
%    useseries       ->   String indicating the type of time series
%                         to use. Can be 'displacement',
%                         'velocity' or 'acceleration'.
  
% Based on function detectsteps.

% Kjartan Halvorsen
% 2003-08-19

[nfr,cols]=size(md);

% Make sure T is even
if mod(T,2)
  T=T-1;
end
T2=T/2;

% Differentiate the data twice
dt = 1/fs;
mdplus=cat(1,md(2:end,:),md(end,:));
mdmin=cat(1,md(1,:),md(1:end-1,:));
mdvel=(mdplus-mdmin)/dt;
mdvel(2:end-1) = mdvel(2:end-1)/2; 

mdvelplus=cat(1,mdvel(2:end,:),mdvel(end,:));
mdvelmin=cat(1,mdvel(1,:),mdvel(1:end-1,:));
mdacc=(mdvelplus-mdvelmin)/dt;
mdacc(2:end-1) = mdacc(2:end-1)/2; 

if strcmp(useseries, 'displacement')
  timeseries = md;
  thesign = 1;
elseif strcmp(useseries, 'velocity')
  timeseries = mdvel;
  thesign = 1;
elseif strcmp(useseries, 'acceleration')
  timeseries = mdacc;
  thesign = -1;
end

% Use the peridogram to estimate the step frequency.
mdfft=fft(md-mean(md),48*nfr);
nperiods=floor(nfr/T); % The approximate number of periods in the data.
psd=abs(mdfft(1:min(100*nperiods,24*nfr)));

maxmdf=max(psd(floor(T/12):end));
peak=find(psd==maxmdf(1));

T=ceil(48*nfr/(peak(1)-1));

% The dominating frequency is found above. Use this to generate a
% sine wave of the same length as the signal. Multiply
% with the signal and take the mean. The result, vv will be a
% function of the phase of the signal.  If we define the
% phase to be zero for a cosine formed signal, then the
% phase will be ph = asin(vv*2);


sw=sin(linspace(0,2*pi*(peak(1)-1)/48,length(timeseries)));

prod=2*mean(sw.*timeseries');
if prod>1
  prod=1;
elseif prod<-1
  prod=-1;
end

ph=fix(asin(prod)/(2*pi)*T);

peaks=(1:T:nfr)'+ph;

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

% Now the peaks are found, find the minimum (maximum) for each period.
np=length(peaks);
stpfrc=zeros(np-1,1);
   
mdc=thesign*timeseries;

mdcfst=mdc((max(1,peaks(1))):(peaks(2)));
mdclst=mdc((peaks(np-1)):(min(nfr,peaks(np))));
mdc = mdc((peaks(2)+1):peaks(np-1));
%keyboard
mdc = reshape(mdc,T,np-3);
   
ffst=find(mdcfst==min(mdcfst));
stpfrc(1)=peaks(1)+ffst(1);

flst=find(mdclst==min(mdclst));
stpfrc(np-1)=peaks(2)+(np-3)*T+flst(1);

   
D=(mdc==kron(ones(T,1),min(mdc)));
% Find the index of the first 1 in each column.
for p=2:np-2
  fi=find(D(:,p-1));
  stpfrc(p)=peaks(2)+(p-2)*T+fi(1);
end
      

%stpfr=stpfrc-1;
stpfr=stpfrc;

