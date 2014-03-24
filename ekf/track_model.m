function [st,dataframes]=track_model(varargin)
% function [st,dataframes]=track_model(bm, mdata, nmrks, filterparam)

% Tracks the model using an ekf.
%
% Input
%    bm        ->   model struct, containing fields 'tws', 'p0'
%                   and 'gcnames'.
%    mdata     ->   Cell array with marker data {attr, md} 
%    nmrks     ->   vector giving the number of markers available 
%		    in each frame.
%    filterparam -> Optional. Acceleration factor. User will be
%                   prompted if this argument is missing.
% Output
%    st        <-   a matrix where each column is a state vector.
%    dataframes<-   a binary vector of length equal to the number of
%                   frames. A 1 indicates that data exists.
%

% Kjartan Halvorsen
% 2003-11-04
%

bm=varargin{1};
mdata=varargin{2};

mnames = getvalue(mdata{1}, 'MARKER_NAMES');

if (nargin > 2 & ~isempty(varargin{3})) 
   nmrks=varargin{3};
else 
   nmrks=length(mnames)*ones(size(mdata{2},1)/3,1);
end

% Prepare marker data.
[initnames, p0, p0vec] = prepare_mdata(bm.p0);
%size(initnames)
%size(p0vec)

%mnames=switchalias(mnames,alilist);
y_observations=extractmarkers(mdata{2}, mnames, initnames);

%TEMPORARY CODE!
y_observations = y_observations(1:80,:);


nfr=size(y_observations,1);

% Replace NaNs with zeros.
y_observations(find(isnan(y_observations)))=0;

% Remove frames at the beginning and end that only contains zeros
%[attr,y_observations,ind_removed,nmrks]=removeemptyframes({},y_observations);
y_observations=y_observations';
dataframes=ones(nfr,1);
%dataframes(ind_removed)=0;

% Enter keyboard mode for debugging
%keyboard

nfr=size(y_observations,2);

% The sampling time
freq=getvalue(mdata{1},'FREQUENCY');
if isstr(freq) freq=str2num(freq); end
if (isempty(freq))
  dt=1/240;
else
  dt=1/freq;
end

% Get filter parameters
% Ask the user to  choose the bandwidth of the tracking filter
if (nargin < 4)
  rrr = input('Enter filter parameter r: ');
  ra=1e3*rrr;
else
  ra=1e3*varargin{4};
end

R2=kron(diag(ones(length(initnames),1)), eye(3));
    
% The form of the process noise covariance is taken from Farina and 
% Studer "Radar data processing". It is the limit as dt*alpha goes
% to zero, or in plain words when the sampling time is much smaller
% than the time constant of the acceleration process.
nst=size(bm.gcnames,1);
ranges = values(bm.gcnames);
Ra=ra*diag(cat(1, ranges{:}));
R1=kron([dt^3/3 dt^2/2; dt^2/2 dt],Ra);

P0=ra*kron(diag([1;10/dt]),eye(nst));

%keyboard


% To find initial state, first convert first frame data to tree
% structure, then call the firststate function
initp = vect2tree(y_observations(:,1), p0);

%x0 = firststate(bm.twists, p0, initp);
x0 = zeros(nst,1);
%keyboard

x0 = cat(1, x0, zeros(size(x0)));
% DEBUG
[y_test, Htest, slsk, y_n] = ...
    observe_mechanism_H(x0, y_observations(:,1), bm.twists, p0, {P0, ...
		    R2});

rmse = y_test - y_observations(:,1);
rmse = sqrt(mean(rmse.*rmse));

%x0=zeros(nst*2,1);


%y_observations(:,1:10)
%p0

disp(['Tracking ', int2str(nst), ' degrees of freedom model ', ...
      'using trajectories of ', int2str(length(initnames)), ' markers...'])

[xhat,Phat]=ekf_mechanism_fixint(y_observations,x0,P0,R1,R2,...
				 dt,bm.twists,p0, rmse);
%[xhat,Phat]=ekf_mechanism_fixlag(y_observations,x0,P0,R1,R2,...
%			          dt,km.twists,p0,6);
   
st=xhat(1:nst,:);

