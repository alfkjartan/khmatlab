function hh=putvalue(varargin)
% h=putvalue(h,key,val) or
% h=putvalue(key,val) (Constructs a new hash)
% Sets (or adds) the value associated with the given key for the hash h.
% Since matlab does not have a hash data structure, the hash is 
% implemented as a nx2 cell array with key and value pairs.

% Kjartan Halvorsen
% 2000-09-12
%
% Revisions
% 2000-11-24   Added code in the case nargin==2: return new hash

if (nargin==2)
   hh=cell(1,2);
   hh{1}=varargin{1};
   hh{2}=varargin{2};
   return
end

h=varargin{1};
key=varargin{2};
val=varargin{3};

hh=h;

if (isempty(hh))
   hh=cell(1,2);
   hh{1}=key;
   hh{2}=val;
   return
end

[p,to]=size(h);

ind=0;
flag=0;
while (~flag & ind<p)
  ind=ind+1;
  if (strcmp(h{ind,1},key))
    hh{ind,2}=val;
    flag=1;
  end
end

if (flag==0)
  hh{p+1,1}=key;
  hh{p+1,2}=val; 
end
