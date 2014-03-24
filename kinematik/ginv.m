function gi=ginv(g)
% function gi=ginv(g)
% Returns the inverse of the rigid body transformation g.

% Kjartan Halvorsen
% 2000-02-27

% Revisions
% 2001-08-05   Handles several frames

[four,four,nfr]=size(g);

if (nfr==1)
   gi=eye(4);
   Ri=g(1:3,1:3)';
   gi(1:3,1:3)=Ri;
   gi(1:3,4)=-Ri*g(1:3,4);
else
   gi=zeros(size(g));

   for i=1:nfr
      Ri=g(1:3,1:3,i)';
      gi(1:3,1:3,i)=Ri;
      gi(1:3,4,i)=-Ri*g(1:3,4,i);
      gi(4,4,i)=1;
   end
end

