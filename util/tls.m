function [X, resnorm, SE, Xlow, Xhigh] = tls(A,B)
% [X, resnorm, SE = tls(A,B)
% Total least squares solution to A*X = B

n = size(A,2); % n is the width of A (A is m by n)
C = [A B]; % C is A augmented with B.
[U S V] = svd(C,0); % find the SVD of C.
%[U S V] = svd(C); % find the SVD of C.
VAB = V(1:n,1+n:end); % Take the block of V consisting of the first n rows and the n+1 to last column
VBB = V(1+n:end,1+n:end); % Take the bottom-right block of V.
X = -VAB/VBB;

if nargout > 1
  % Return also the frobenius norm of the residuals (errors)
  VABBB = V(:,1:n:end);
  EF = - C*VABBB*VABBB';
  resnorm = norm(EF, 'fro');
end

if ( (nargout > 2) ) % & (length(X) == 1) )
  % Return the standard error of the first (slope) estimate
  
  % Residuals
  E = B - A*X;
  n = size(A,1);
  SE = sqrt(sum(E.^2, 1) / (n-2)) ...
       / sqrt( sum((A(:,1) - mean(A(:,1))).^2 ,1) );
  

  % The high value of the conf interval
  Xhigh = X + 2*SE;
  
  % The low value of the conf interval
  Xlow = X - 2*SE;
  

end