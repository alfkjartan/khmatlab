function tw=expcoord(g)
% function tw=expcoord(g)
% Computes the twist (exponential coordinates) associated with the
% (4 x 4) rigid transformation g.

% Kjartan Halvorsen
% 1999-05-31

thr=1e12;

R=g(1:3,1:3);
p=g(1:3,4);
niv=natinvar(R);
w=niv(1:3);
theta=niv(4);
w_hat=vect2mat(w);
A=(eye(3)-R)*w_hat+w*w'*theta;
if (cond(A)>thr)
  theta=norm(p);
  v=p/norm(p);
else
  v=inv(A)*p;
end
tw=twist([v;w]);
tw.theta=theta;
tw=twist(tw);