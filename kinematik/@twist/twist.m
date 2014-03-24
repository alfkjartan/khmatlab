function tw=twist(varargin)
% Constructor for twist objects.
% Usage:
% tw=twist([v;omega],theta)
% tw=twist([R d;0 0 0 1])
% tw=twist(tw)
%
% Reference: Murray, Li, Sastry: A mathematical introduction to
% robotics. 

% Kjartan Halvorsen
% 1999-05-31

tolerance=1e-8;

if (nargin==0)
  tw.coordinates=zeros(6,1);
  tw.theta=0;
else
  
  if (isa(varargin{1},'twist'))
    tw=varargin{1};
    return
  else
    vw=varargin{1};
    vw=vw(:);
    if (length(vw(:))==16) % assume rigid motion is given
      tw=expcoord(varargin{1});
      return
    else
      v=vw(1:3);
      w=vw(4:6);
    end
    if (nargin==1)
      th=0;
    else
      th=varargin{2};
    end
  end


  if(norm(w)~=0)
    if (abs(norm(w)-1)>tolerance)
      warning('the omega vector should have unit length')
      w=w/norm(w);
    end
  end

  tw.coordinates=[v;w];
  tw.theta=th;
end

tw=class(tw,'twist');