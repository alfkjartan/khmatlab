function [np,okfrs] = movePoints(g, p)
%  np = movePoints(g, p)
% Returns the paths of points in p, transformed according to g
%
% Input
%   g          ->   (4 x 4 x nfr) or (nfr x 16) sequence of rigid
%                   body transformations.
%   p          ->   Set of 3d points. (4 x n) or (3 x n).
% Output
%   np         <-   Paths of 3d points. (3 x n x nfr)

% Kjartan Halvorsen
% 2003-11-17

% Revisions
% 2011-01-17
%              Will not rotate marker positions that are (0,0,0)

if (size(p, 1) == 3)
  p = cat(1, p, ones(1, size(p, 2)));
end

if (size(g, 2) == 16) 
  nfr = size(g, 1);
else
  nfr = size(g, 3);
end

np = zeros(3, size(p, 2), nfr);

allfrs = 1:nfr;
okfrs = find (sum(p(1:3,:)) ~= 0);

if (size(g, 2) == 16) 
  for f=allfrs
    npp = zeros(size(p));
    gg = reshape(g(f,:), 4, 4);
    npp(:,okfrs) = gg*p(:,okfrs);
    np(:,:,f) = npp(1:3, :);
  end
else
  for f=allfrs
    npp = zeros(size(p));
    npp(:,okfrs) = g(:,:,f)*p(:,okfrs);
    np(:,:,f) = npp(1:3, :);
  end
end
