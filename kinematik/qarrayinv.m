function qi = qarrayinv(q)
% Computes the inverse quaternion of each quernion (rows) of q

% Kjartan halvorsen
% 2012-10-24

qi = zeros(size(q));

for i=1:size(q,1)
  qi(i,:) = qinv(quaternion(q(i,1), q(i,2), q(i,3), q(i,4)));
end
