function M=vect2mat(v)
% function M=vect2mat(v)
% Returns the skew-symmetric matrix associated with v.

% Kjartan Halvorsen
% 1999-05-31

v_hat=zeros(3,3);
v_hat(1,2)=-v(3);
v_hat(1,3)=v(2);
v_hat(2,1)=v(3);
v_hat(2,3)=-v(1);
v_hat(3,1)=-v(2);
v_hat(3,2)=v(1);

M=v_hat;