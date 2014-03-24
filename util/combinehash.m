function hh=combinehash(h1,h2)
%  hh=combinehash(h1,h2)
% Combines the two hashes ((n x 2) cell arrays) h1 and h2.
% elements with equal keys will not be duplicated, they
% will come from h1.

% Kjartan Halvorsen
% 2002-03-24

hh=h1;

hlength=size(h1,1);
h1keys=h1(:,1);

for l=1:length(h2(:,1))
   slask=setdiff(h2{l,1},h1keys);
   if ~isempty(slask)
      hh(hlength,:)=h2(l,:);
      hlength=hlength+1;
   end
end
