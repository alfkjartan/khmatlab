function [th, d] = invkin_singleaxis(p, q, w)
%  [th, d] = invkin_singleaxis(p, q, w)
% 
% Solves 
%     q'= Rp' = exp{hat(w)*th}p',
% where q' is q with the centroid subtracted.
% d = q' - Rp'
%
% input
%    p     ->  (3 x m) set of markers
%    q     ->  (3 x m) set of markers
%    w     ->  (3 x 1) unit vector in the direction of the axis.
  
% Kjartan Halvorsen
% 2004-03-16

G = soder(cat(1, (p(:))', (q(:))'));

th = w'*hat(G(1:3,1:3));

d = G(1:3,4);



  
 

  