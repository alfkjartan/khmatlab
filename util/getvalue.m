function v=getvalue(h,key)
% Returns the value associated with the given key for the hash h.
% Since matlab does not have a hash data structure, the hash is 
% implemented as a nx2 cell array with key and value pairs.

% Kjartan Halvorsen
% 2000-09-12

%%DEBUG
%key=key
  
% Set to 1 if function should warn when element not found
warna = 0;
  
[p,to]=size(h);
v='';
ind=0;
flag=1;
while (isempty(v) & ind<p)
  ind=ind+1;
  if (strcmp(h{ind,1},key))
    v=h{ind,2};
    flag=0;
  end
end

if (warna & flag)
   warning(['Key ',key,' not found'])
end
