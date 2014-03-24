function [np,res]=getRelMotion(distp,proxp)
% function np=getRelMotion(distp,proxp)
% returns the motion of distp as seen from proxp.
%
% Input
%    distp    ->    Data for distal markers. (nfr x 3*ndmarks)
%                   matrix
%    proxp    ->    Data for proximal markers. (nfr x 3*npmarks)
%                   matrix
% Output
%    np       <-    Data for distal markers as observed in a
%                   coordinate system fixed in the proximal segment
%                   (nfr x 3ndmarks) matrix
  
% Kjartan Halvorsen
  
[T,res]=getMotion(proxp);
[m,n]=size(T);
[md,nd]=size(distp);
[mp,np]=size(proxp);
labels=mod(nd,3);
points=floor(nd/3);

if (mp~=md)
   error('data sets distp and proxp must have the same length')
end

dp=distp;
if (labels==1)
   dp(:,1)=[];
end

np=zeros(m,points*3);

ppp = ones(4, points);
for i=1:md
   pp=dp(i,:);
   if (~hasmissing(pp) & ~hasmissing(proxp(i,:))) % no missing markers
     ppp(1:3,:)=reshape(pp,3,points);
     TT=reshape(T(i,(labels+1):(16+labels)),4,4);
     TTT=inv(TT);
     pppp=TTT*ppp;
     pppp(4,:)=[];
     np(i,1:(points*3))=reshape(pppp,1,points*3);
   end
end

if (labels==1)
	np=[distp(:,1) np];
end
