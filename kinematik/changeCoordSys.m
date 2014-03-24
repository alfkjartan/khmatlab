function np = changeCoordSys(g, p)
%  np = changeCoordSys(g, p)
% Returns the paths of points in p, transformed to a coordinate
% system moving as described by g
%
% Input
%   g          ->   (4 x 4 x nfr) or (nfr x 16) sequence of rigid
%                   body transformations.
%   p          ->   Paths of 3d points. (3 x n x nfr) 
% Output
%   np         <-   Paths of 3d points. (3 x n x nfr)

% Kjartan Halvorsen
% 2003-11-17
  
if (size(g, 3) == 1) 
  nfr = size(g, 1);
else
  nfr = size(g, 3);
end

if (size(p,3) ~= nfr)
  error('Wrong number of frames in the two input arguments')
end

nps = size(p,2);

np = zeros(3, size(p, 2), nfr);

if (size(g, 3) == 1) 
  for f=1:nfr
    gg = reshape(g(f,:), 4, 4);
    npp = ginv(gg)*cat(1, p(:,:,f), ones(1, nps)) ;
    np(:,:,f) = npp(1:3, :);
  end
else
  for f=1:nfr
    npp = ginv(g(:,:,f)) * cat(1,  p(:,:,f), ones(1, nps)) ;
    np(:,:,f) = npp(1:3, :);
  end
end
  