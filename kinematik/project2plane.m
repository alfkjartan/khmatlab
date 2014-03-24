function nmd = project2plane(md, plane_normal, vertical)
% nmd = project2plane(md, plane_normal, vertical)
% Projects the marker data onto a plane. 
%
% Input
%    md           ->   marker data (n x 3m).
%    plane_normal ->   unit vector normal to 2d plane
%    vertical     ->   unit vector in the vertical 2d
%                      direction. That is vertical in the 2d view.
% Output
%    nmd          <-   Transformed marker data. (n x 3m) homogenous
%                      2d marker data.
  
% Kjartan Halvorsen
% 2004-03-23

if (nargin > 0)
  e_z = plane_normal;
  e_y = vertical;
  e_y = e_y - (e_y'*e_z)*e_z;
  e_y = e_y / norm(e_y);

  % Rotation matrix that rotates world coordinates into 2d plane coordinates
  R = (cat(2, cross(e_y, e_z), e_y, e_z))';

  [n,m3] = size(md);
  nmd = R * reshape(md', 3, m3/3*n);
  nmd(3,:) = 1;
  nmd = reshape(nmd, m3, n);
  nmd = nmd';

else % Unit test
end

  % Set up test case
  
  


