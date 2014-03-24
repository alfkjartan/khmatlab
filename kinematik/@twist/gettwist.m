function t=gettwist(tw)
% function t=gettwist(tw)
% Returns the (4 x 4) se(3) element twist.

% Kjartan Halvorsen

vw=tw.coordinates;
v=vw(1:3);
v(:);
w=vw(4:6);

w_hat=vect2mat(w);
t=[w_hat v;0 0 0 0];