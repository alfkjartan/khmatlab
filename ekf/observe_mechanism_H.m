function [y,H,Gs,ytn]=observe_mechanism_H(x,yt,tws,p0,Gprox,dgdx)
%  [y, H, Gs, ytn]=observe_mechanism_H(x, yt, tws, p0s, {Pt, R2})
% "Observation function" for a general mechanism
% (thigh and shank e.g.).
% Given the state of the system, the twists and the position of the markers
% in the reference (first) frame it returns the current position of
% the same markers.
%
% Input
%    x     ->   the current state. (2n x N) vector. 
%    yt    ->   the measurements.(3*m x 1) vector. Used for
%               detecting missing markers.
%    tws   ->   nested cell array with twists
%    p0s   ->   nested cell array with markers (2 x 1)
%    {Pt,R2}  -> cell array containing the estimation error
%                covariance matrix and the measurement noise
%                covariance matrix. Optional.
% Output
%    y     <-   the current marker positions. (3*m x 1) vector.
%    H     <-   Jacobian or the linearized observation matrix.
%               (3m x 2n) matrix. 
%    Gs    <-   Transformation of the root segment. Normally only
%               used in the recursive calls.
%    ytn   <-   if Pt and R2 are provided, contains the corrected observation
%               vector with some elements swopped.

% Kjartan Halvorsen
% 2002-05-30
%
%
% Revisions
% 2004-10-26   Added a check for switched identitity of markers and
%              for outliers.
% 2009-06-24   Changed to use NaN to identify missing markers. This
%              is in order to be consistent with soder

% Gprox contains the rigid body motion of the proximal segment.
% dgdx is a 3 dimensional array containing the derivative of the
% proximal rigid body transformation wrt the states.

% gsprox is a 3 dimensional array containing the exponential of all
% the proximal twists.
% twsprox is a cell array with the twists defining the motion of
% the proximal segment. xprox contains the state (generalized
% coordinates) for the same twists.
% This function is called recursively to calculate y and H

%%try
  
if nargin==4  % Initial call
%  n=length(x);
%  x=x(1:n/2,:);
  Gprox=eye(4);
  dgdx=[];
elseif nargin==5  % Initial call, with check for errors in marker data
  Pt = Gprox{1};
  R2 = Gprox{2};
  Gprox = eye(4)
  dgdx = [];
  
  % Confidence level for the size of innovations. Taken from a
  % table of the Xi2 distribution with 3 degrees of freedom. 
  % P(Z < 7) = 0.95
  dlim = 7;
  
  % P(Z < 1) = 0.2
  dlim2 = 1;
end

% First the "own" marker positions and Jacobian, then the branches

mytws=tws{1};
nn=length(mytws);


myx=x(1:nn);

if isempty(dgdx)
  nprox=0;
else
  nprox=size(dgdx,3);
end

myp0=p0{1};
[lep,mp]=size(myp0);
if lep==3
  myp0=cat(1,myp0,ones(1,mp));
end 

gs=zeros(4,4,nn);
xsis=zeros(4,4,nn);

gg=eye(4);
for st=1:nn
  %gs(:,:,st)=expr(mytws{st},x(st));
  gs(:,:,st) = expm(mytws{st} * x(st));
  xsis(:,:,st)=mytws{st};
  gg=gg*gs(:,:,st);
end

% The position of the markers
if ~isempty(myp0)
  pp=Gprox*gg*myp0;
  pp(4,:)=[];
  y=pp(:);
else
  y=[];
end

if nargout>2
  Gs={Gprox*gg};
end

%return
% ------------------------------------------------
%  The Jacobian
% ------------------------------------------------


% transformations upto angle and after angle
if (nn>0)
  before=zeros(4,4,nn);
  after=zeros(4,4,nn);

  before(:,:,1)=Gprox*gs(:,:,1);
  after(:,:,end)=eye(4);
  for k=1:nn-1
    before(:,:,k+1)=before(:,:,k)*gs(:,:,k+1);
    after(:,:,end-k)=gs(:,:,end-k+1)*after(:,:,end-k+1);
  end
end

mydgdx=zeros(4,4,nprox+nn);
for k=1:nprox
  mydgdx(:,:,k)=dgdx(:,:,k)*gg;
end
for k=nprox+1:nprox+nn
  mydgdx(:,:,k)=before(:,:,k-nprox)*xsis(:,:,k-nprox)*after(:,:,k-nprox);
end


%mydgdx
%before
%after
%xsis

%pause

if ~isempty(myp0)
  dpdphi_p=zeros(4,(nprox+nn)*mp);
  for k=1:(nn+nprox)
    dpdphi_p(:,(k-1)*mp+1:k*mp)=mydgdx(:,:,k)*myp0;
  end
  dydphi_p=[eye(3) zeros(3,1)]*dpdphi_p;

  HH{1}=reshape(dydphi_p,mp*3,nprox+nn);
else
  HH{1}=[];
end


% ----------------------------------------------------------
% Before finishing the job, call function recursively
% to get the position and jacobian of the distal segments.
% ----------------------------------------------------------

if (length(tws)>1) % Branches exist
%% For debug  warning(['Found ', int2str(length(tws)-1), ' branches'])
  twsbr=tws(2:length(tws));
  ndist=0;
  if (isempty(yt))
    ytdist = [];
  else
    ytdist=yt( (length(y)+1):length(yt),1);
  end
  xelems = size(x,1);
  for br=1:length(twsbr)
    [ydist,Hdist,Gdists]=observe_mechanism_H(x((nn+ndist+1):xelems, 1),...
				      ytdist,...
				      twsbr{br},...
				      p0{br+1},...
				      Gprox*gg,...
				      mydgdx);
    y=[y;ydist];

    if nargout==3
      Gs=cat(1,Gs,Gdists);
    end
    
    if isempty(Hdist)
      HH{br+1}=[];
    else
%     size(ydist)
%     size(Hdist)
%     br
%     p0{br+1}
    
      Hdistcols = size(Hdist, 2);
      HH{br+1}=[Hdist(:,1:nprox+nn) zeros(length(ydist),ndist) ...
	  	Hdist(:,nprox+nn+1:Hdistcols)];
      
      ndistbr=size(HH{br+1},2)-size(HH{1},2);
      for bl=1:br
	      HH{bl}=[HH{bl} zeros(size(HH{bl},1),ndistbr)];
      end  
      ndist=ndist+size(Hdist,2)-nn-nprox;
    end  
  end
  
  %% Construct the Jacobian
  H=cat(1,HH{:});

else
  H=HH{1};
end

if nargin<6 % Initial call
   H=[H zeros(size(H))];
   ind=find(yt==0);
   %%ind=find(isnan(yt));
   H(ind,:)=0;
   
   if (nargin==5)
     %% Check for markers witch may be erronous (innovation too
     %% large), and markers which may have swithed identity.
     
     %% The innovations
     ytilde = y - yt
     
     % The covariance matrix of the innovations
     Q = H*Pt*H' + R2;
     
     nmrks = length(yt)/3;
     suspect = [];
     for m=1:nmrks
       d = ytilde((m-1)*3+1:m*3)' * ...
	         inv(Q((m-1)*3+1:m*3, (m-1)*3+1:m*3)) * ...
	         ytilde((m-1)*3+1:m*3)

       if d > dlim
	       suspect = cat(1, suspect, m);
       end 
     end 
     
     if ~isempty(suspect)
       warning(['Found markers ', int2str(suspect'), ...
		' to be erronous']);
       
       % Compute a matrix containing the Mahalanobis distance
       % between each of the suspected predicted markers and
       % observations.
       dmatrix = repmat(1e10,length(suspect));
       for i = 1:length(suspect)
	       m=suspect(i);
	       Qminv = inv(Q((m-1)*3+1:m*3, (m-1)*3+1:m*3));
	       for j = 1:length(suspect)
	         if (j ~= i)
	           k=suspect(i);
	           dmatrix(i,j) = (y((m-1)*3+1:m*3) - yt((k-1)*3+1:k*3))' ...
		             * Qminv * ...
	               (y((m-1)*3+1:m*3) - yt((k-1)*3+1:k*3));
	         end 
	       end 
       end 

       %% If symmetric elements of dmatrix have small values, then
       % the markers have been swopped. 
       % Try to find best matches first.

       dmins = min(dmatrix,[], 2);
       [slsk, minind] = sort(dmins);
       
       taken = [];
       swoptable = [];
       for i = 1:length(minind)
	       ind = minind(i);
	       if dmins(ind) < dlim2 % Likely point match
	         colind = find(dmatrix(ind,:) == dmins(ind));
	         if ~ismember(colind, taken)
	           taken = cat(1, taken, colind);
	           swoptable = cat(1, swoptable, ...
			                       [suspect(ind) suspect(colind)]);
	         end
	       end
       end
       
       % Now swop the data vector yt and return in ytn.
       ytn = yt;
       for i=1:size(swoptable, 1)
	       k = swoptable(i,1);
	       l = swoptable(i,2);
	       ytn((k-1)*3+1:k*3) = yt((l-1)*3+1:l*3); 
       end
     
       % For the rest of the suspected markers that have not been
       % corrected, set corresponding rows of H to zero
       if isempty(swoptable)
	       notcorrected = suspect;
       else
	       notcorrected = setdiff(suspect, swoptable(:,1));
       end 
       
       for nc = notcorrected'
	       H((nc-1)*3+1:nc*3,:) = 0;
       end 
     else
       ytn = yt;
     end 
   end 
 end
 
%%catch
%%  keyboard
%%end
 
