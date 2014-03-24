function vs = rotate(vb, Rs, Rsb) 
%% vs = rotate(vb, Rs, Rsb) 
%% Returns the vectors in v(i,:) rotated by R*R0, that is
%%   vn(i,:) = v(i,:) * Rsb'*Rs(i)

%% Kjartan Halvorsen
%% 2012-05-09

vs = zeros(size(vb));

nfr = min(size(vs,1), size(Rs,1));

Rbs = Rsb';

for i=1:nfr
  Rsi = reshape(Rs(i,1:12), 4, 3);
  vs(i,:) = vb(i,:)*Rbs*Rsi(1:3,1:3)';
end
