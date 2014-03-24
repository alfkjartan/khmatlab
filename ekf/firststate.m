function st=firststate(tws, p0, firstp, g_prox)
% function st=firststate(tws, p0, firstp, g_prox)
% Finds generalized coordinates that approximate the position of
% markers  state for the kinmodel that matches the given position of 
% the joint centers.
%
% Input
%    tws      ->   Twists, tree structure
%    p0       ->   (tree) 3d coordinates of markers at reference pos
%    firstp    ->   (tree) 3d coordinates of markers at first frame
%    g_prox   ->   pose of proximal object.
% 
% Output
%    st       <-   states of twists
%

% Kjartan Halvorsen
% 1999-09-20

  
if nargin==3 % Initial call
  g_prox=eye(4);
end

invg=ginv(g_prox);

mytws=tws{1};
nn=length(mytws);

  
myp0=p0{1};
[lep,mp]=size(myp0);
if lep==3
  myp0=cat(1,myp0,ones(1,mp));
end
myp = myp0;

myfirstp = firstp{1};

[lep,mp]=size(myfirstp);
if lep==3
  myfirstp=cat(1,myfirstp,ones(1,mp));
end

st = [];

if isempty(myfirstp)
    warning(['cannot solve inverse kinematics problem. ',...
	     'No markers.'])
    st=zeros(length(mytws),1);
else
  myfirstp = invg*myfirstp;
  if (mp > 2)
    % Compute transformation, call g2var
    myg = soder(cat(1, reshape(myp(1:3,:), 1, mp*3),...
    		    reshape(myfirstp(1:3,:), 1, mp*3)));
    
    if ~isnan(sum(sum(myg)))
      st = g2var(mytws, myg);
    else
      warning(['cannot solve inverse kinematics problem. ',...
	       'Transformation is NaN.'])
      st = zeros(length(mytws),1);
    end
  
    if (nn==3)
     % keyboard
    end
    
  else
    if (mp ==2)
      %two points. Use the middlepoint if not one is NaN
      if ~isnan(sum(sum(myfirstp)))
	mypm = mean(myp(1:3,:),2);
	myfirstpm = mean(myfirstp(1:3,:),2);
      else
	if isnan(sum(myfirstp(:,1)))
	  if isnan(sum(myfirstp(:,2)))
n	    % OBS. No data available.
	    st = zeros(length(mytws),1);
	  else
	    mypm =myp(1:3,2);
	    myfirstpm = myfirstp(1:3,2);
	  end
	else % The other point must be NaN, but the first ok
	  mypm =myp(1:3,1);
	  myfirstpm = myfirstp(1:3,1);
	end
      end % if ~isnan(...)
    else % Only one point. Check if nan
      if isnan(sum(myfirstp(:,1)))
	% OBS. No data available.
	st = zeros(length(mytws),1);
      else
	mypm =myp(1:3,1);
	myfirstpm = myfirstp(1:3,1);
      end
    end
    
    % Single point or two points only, in this case there should
    % either be one or two twists. 
  
    if isempty(st)
      if (nn == 1)
	[st,fl]=padkah1(mytws{1},mypm, myfirstpm);
      elseif (nn == 2)
	[st,fl]=padkah2(mytws{:},mypm,myfirstpm);
      else
	warning(['cannot solve inverse kinematics problem. ',...
		 'Setting angles to zero.'])
	st=zeros(length(mytws),1);
      end
    end
  end
end

if length(tws)>1 % Branches exist

  try
  gg=g_prox;
  for i=1:nn
   gg = gg * expm(mytws{i}*st(i));
  end
  
  if (nn==6)
    %%disp('firststate')
    %%myg
    %%gg
    %keyboard
  end
  catch
    keyboard
  end
  
  twsbr=tws(2:end);
  for br=1:length(twsbr)
    [stdist]=firststate(twsbr{br}, p0{br+1}, firstp{br+1},...
			gg);
    st = cat(1, st, stdist);
  end
end
