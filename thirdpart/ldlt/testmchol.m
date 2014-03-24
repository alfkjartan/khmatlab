%
% m-file to test the ldlt m-files.
%
A=rand(5,5);
while (rank(A) ~= 5),
  A=rand(5,5);
end;
A=A*A';
%
%  Testing ldlt.
%
[L,D]=ldlt(A);
disp('ldlt test: should be close to 0')
norm(A-L*D*L')
%
%  testing ldltup.
%
v=rand(5,1);
[newL,newD]=ldltup(L,D,v);
disp('ldltup test: should be close to 0')
norm(A+v*v'-newL*newD*newL')
%
% ldltdown test.
%
w=v/1000;
[newL,newD]=ldltdown(L,D,w);
disp('ldltdown test: should be close to 0')
norm(A-w*w'-newL*newD*newL')
%
% mchol test.
%
[L,D,E,pneg]=mchol(A);
disp('mchol test 1: should be close to 0')
norm(A-L*D*L')
B=A-min(eig(A))*1.5*eye(5);
[L,D,E,pneg]=mchol(B);
disp('mchol test 2: should be close to 0')
norm(B+E-L*D*L')
disp('mchol test 3: should be negative')
pneg'*B*pneg
%
%
%
%
