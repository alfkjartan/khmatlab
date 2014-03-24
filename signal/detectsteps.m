function [stpfr,mdvelw,isreflex]=detectsteps(md,T)
% stpfr=detectsteps(md,T)
% Tries to identify steps in the marker data. md contains only vertical 
% coordinates of markers that may give information about the steps. 
% The coordinate data are filtered, then differentiated. The start
% of the step is taken to be at the lower peak of the velocity,
% i.e. when the downward velocity changes dramatically. This
% implies the positive direction being upwards.
%
% Based on earlier function with the same name.
%
% Input
%    md              ->   marker data (N x m)
%    T               ->   Approximate period length

% Kjartan Halvorsen
% 2002-12-12

reflexthr=0.3;
isreflex = 0;

[nfr,cols]=size(md);

% Make sure T is even
if mod(T,2)
  T=T-1;
end
T2=T/2;

smplfreq=nfr;
cutofffreq=1.5*nfr/T;

% Filter the data hard
mdd=detrend(md);
mdf=filtermdata(mdd,smplfreq,cutofffreq,4);

% And less hard
mdfw=filtermdata(mdd,smplfreq,smplfreq/24,4);

% Differentiate the data
mdplus=cat(1,mdf(2:end,:),mdf(end,:));
mdmin=cat(1,mdf(1,:),mdf(1:end-1,:));
mdvel=mdplus-mdmin;

mdplus=cat(1,mdfw(2:end,:),mdfw(end,:));
mdmin=cat(1,mdfw(1,:),mdfw(1:end-1,:));
mdvelw=mdplus-mdmin;

% Check if velocity below threshold, if so, reflex file
stvel=std(mdvelw);
if (stvel<reflexthr)
  isreflex = 1;
  stpfr=[];
  return
elseif (stvel<10*reflexthr)
  answer=questdlg('Is this a reflex or stance file?',...
		  'Question',...
		  'Yes','No','Yes');
  if strcmp(answer,'Yes')
    isreflex = 1;
    stpfr=[];
    return
  end
end

% Use the peridogram to estimate the step frequency.
mdfft=fft(mdd,48*nfr);
nperiods=floor(nfr/T); % The approximate number of periods in the data.
psd=abs(mdfft(1:min(100*nperiods,24*nfr)));

maxmdf=max(psd);
peak=find(psd==maxmdf);

T=ceil(48*nfr/(peak(1)-1));

% Do the following: 
% The dominating frequency is found above. Use this to generate a
% sine wave of the same length as the velocity signal. Multiply
% with the velocity signal and take the mean. The result, vv will be a
% function of the phase of the velocity signal.  If we define the
% phase to be zero for a cosine formed velocity signal, then the
% phase will be 
% ph = asin(vv*2);

sw=sin(linspace(0,2*pi*(peak(1)-1)/48,length(mdvel)));

prod=2*mean(sw.*mdvel');
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

% Now the peaks are found, compute the minimum for each period.
np=length(peaks);
stpfrc=zeros(np-1,1);
   
mdc=mdvelw;
      
mdcfst=mdc((max(1,peaks(1))):(peaks(2)));
mdclst=mdc((peaks(np-1)):(min(nfr,peaks(np))));
mdc=reshape(mdc((peaks(2)):(peaks(np-1)-1)),T,np-3);
   
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
      

stpfr{1}=stpfrc-1;

%   catch
%      uiwait(warndlg(['Unable to find enough peaks in the data.', ...
%                      'Try to detect the steps manually']));
%      stpfr={};
%      return
%   end
      
   