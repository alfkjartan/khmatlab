function [ib,ie,iz]=beginEnd(p)
%function [ib,ie,iz]=beginEnd(p)
% Returns the index of the first and last row with data. Returns in iz also
% the indices of frames between ib and ie that contain zeros

% Kjartan Halvorsen
% 2000-05-15

[m,n]=size(p);

if (n>1)
  pp=(abs(p)>0);
  tmp1=(sum(pp'))';
  tmp2=(tmp1==n);
else
  tmp2=(abs(p)>0);
end

ib=-1;
ie=-1;
k=0;
while(ib<0)
  k=k+1;
  if (ib<0 & tmp2(k)==1) ib=k; end
end

k=m+1;
while(ie<0)
  k=k-1;
  if (ie<0 & tmp2(k)==1) ie=k; end
end

iz=ib-1+find(~tmp2(ib:ie));
