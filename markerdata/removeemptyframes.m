function [nattr,nmd,frs,nmarks]=removeemptyframes(attr,md)
%  [nattr,nmd,frs,nmarks]=removeemptyframes(attr,md)
% Removes the frames at the beginning and the end which only contain
% zeros. 
%
% Input
%      attr      ->   The attribute hash
%      md        ->   The marker data (nfr x m*3)
% Output
%      nattr     <-   The new attribute hash
%      nmd       <-   The new marker data
%      frs       <-   The index of the removed frames.
%      nmarks    <-   A (nfr x 1) vector giving the number of markers
%		      available at all frames.
%

% Kjartan Halvorsen
% 2000-12-14

% Extract the part of the file that contains data. 

ms=length(md(1,:))/3;

ib=1;
ie=length(md(:,1));

nmarks=zeros(ie,1);

for m=1:ms
   [ibm,iem,iz]=beginEnd(md(:,(m-1)*3+1:m*3));
   ie=min(ie,iem);
   ib=max(ib,ibm);
   nmarks(ibm:iem)=nmarks(ibm:iem)+1;
   nmarks(iz)=nmarks(iz)-1;
end

frs=setdiff((1:length(md(:,1))),(ib:ie));
nmd=md(ib:ie,:);
nattr=putvalue(attr,'NO_OF_FRAMES',length(nmd(:,1)));
