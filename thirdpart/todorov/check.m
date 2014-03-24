%err = check(segment, [parms])
%
% Check observability
%
% segment: kinematic model (see residual.m for details)
% parms: optional matrix, each column is a set of values for the
%        constant (w) parameters of the kinematic model
%
% err: 3x1 cell array with elements corresponding to
%        (estimation, one-sample sysid, full sysid)
%      each element gives the null-space of the corresponding
%      observability matrix; if empty, the model is observable

% Copyright (C) Emanuel Todorov, 2006-2007

function Err = check(segment, varargin)

stdPos = 20;                              % stdev position
stdRot = 1;                               % stdev rotation
nRep = 2;                                 % repetitions
factor = 2;                               % body states = factor * nW
tol = 1E-10;                              % tolerance for rank computation

[map, info] = prepare(segment,3);         % prepare model
nDat = max(10, factor * map.nW);

flagP = (info.type(1:map.nP+map.nW)==1);  % flag p-elements
flagW = (info.type(1:map.nP+map.nW)==2);  % flag w-elements

if ~isempty(varargin),                    % user-supplied model parameters
   nRep = size(varargin{1},2);
   xw = zeros(map.nX,nRep);
   xw(info.type==2,:) = varargin{1};
else
   xw = randn(map.nX,nRep)*stdPos.* ...   % random model parameters
         repmat((info.type==2 & info.spatial==1), [1,nRep]) + ...
        randn(map.nX,nRep)*stdRot.* ...
         repmat((info.type==2 & info.spatial==2), [1,nRep]);
end

szBig = map.nP*nDat+map.nW;               % size of BIG matrix
iEnd = (map.nP*nDat+1:szBig)';            % index of end submatrix
Err = cell(3,1);                          % clear errors
opts = struct('issym',1,'isreal',1,'disp',0);  % eigs options


for rep = 1:nRep                          % loop over parameter sets
   BIG = sparse(szBig,szBig);             % allocate BIG matrix

   for dat = 1:nDat                       % loop over body states
      x = xw(:,rep) + ...                 % generate data point
          randn(map.nX,1)*stdPos.*(info.type==1 & info.spatial==1) + ...
          randn(map.nX,1)*stdRot.*(info.type==1 & info.spatial==2);

      [r,y] = residual(zeros(map.nY,1), x, segment, map);
      [r,y,J] = residual(y, x, segment, map);      % compute Jacobian

      JP = J(:,flagP);                    % Jacobian of p-elements
      JW = J(:,flagW);                    % Jacobian of w-elements

      if rank(JP,tol)<map.nP & isempty(Err{1}),
         Err{1} = null(JP);               % estimation failure
      end

      if rank(JW,tol)<map.nW & isempty(Err{2}),
         Err{2} = null(JW);               % trivial sysid failure
      end

      iSub = (map.nP*(dat-1)+1:map.nP*dat)';  % index of current submatrix
      BIG(iSub,iSub) = JP'*JP;            % fill BIG
      tmp = JP'*JW;
      BIG(iSub,iEnd) = tmp;
      BIG(iEnd,iSub) = tmp';
      BIG(iEnd,iEnd) = BIG(iEnd,iEnd) + JW'*JW;
   end

   if map.nW>0 & isempty(Err{3}),
      BIG = speye(szBig,szBig) - BIG;        % sparse eigensolver
      [V,D] = eigs(BIG,10,'la',opts);
      bad = find(abs(1-diag(D))<tol);
   
      if ~isempty(bad),                      % sysid failure
         V = V(end-map.nW+1:end,bad);
         Err{3} = orth(V);
      end
   end
end


names = {'Estimation', 'One-sample SysID', 'Full SysID'};
for i=1:3                                    % print error info
   if isempty(Err{i}),
      fprintf('%s OK\n', names{i});
   else
      fprintf('%s impossible, null-space is %d-dimensional\n', ...
         names{i}, size(Err{i},2));
   end
end
