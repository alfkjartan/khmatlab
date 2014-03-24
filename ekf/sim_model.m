function [m,mnames, xpoints, xpnames, H]=sim_model(km, states, extrapoints)
%  [m, mnames, xpoints, H]=sim_model(km,states,extrapoints)
% Returns the motion of the segments in the model given the sequence 
% (matrix) of state vectors. 
%
% Input
%    km       ->   The model struct.
%    states   ->   The sequence of states, (nsts x nfr)
%    extrapoints ->   string identifying extra field of km with
%                     points to simulate.
% Output
%    m        <-   Simulated marker data
%    mnames   <-   Cell array with the names of the simulated markers
%    xpoints  <-   simulated extra points. If 'CoM', then the
%                  function will also look for a field 'mass', and 
%                  provide the total CoM as well as the individual
%                  CoM for each segment.
%    xpnames   <-  Cell array with the names of the extra simulated markers
%    H         <-  (3*nmrks x nsts x nfr) The Jacobian matrix.


% Kjartan Halvorsen
% 2003-02-05
%
% Revisions
% 2009-06-25   Added calculation of total CoM if input argument
%              'extrapoints' equals 'CoM'.

[nst,nfr]=size(states);

% Prepare marker data.
[mnames, p0] = prepare_mdata(km.p0);

ys=zeros(3*length(mnames),nfr);
H = zeros(3*length(mnames), nst, nfr);
for i=1:nfr
  [ys(:,i),Hi]=observe_mechanism_H([states(:,i); zeros(nst,1)],[], ...
			      km.twists, p0);
  H(:,:,i) = Hi(:,1:nst);
end
m=ys';

if ( nargin==3 & isfield(km, extrapoints) )
  xp = getfield(km, extrapoints);
  [xpnames, xp] = prepare_mdata(xp);
  ytest=observe_mechanism_H([states(:,1); zeros(nst,1)],[], ...
			      km.twists, xp);
  
  ys=zeros(length(ytest),nfr);
  for i=1:nfr
    [ys(:,i)]=observe_mechanism_H([states(:,i); zeros(nst,1)],[], ...
				  km.twists, xp);
  end

  if strcmp(extrapoints, 'CoM')
    if isfield(km, 'mass')
      mass = tree2vect(km.mass);
      ns = length(mass);
      mass3 = diag(kron(mass, ones(3,1)));
    
      yweighted = mass3*ys;
      ytot = cat(1, ...
		 sum(yweighted(1:3:end,:)),...
		 sum(yweighted(2:3:end,:)),...
		 sum(yweighted(3:3:end,:)));
	       
      xpnametot = {'CoM'};
    else
      ytot=[];
      xpnametot = {};
    end
  else
      ytot=[];
      xpnametot = {};
  end

  
  xpoints=(cat(1, ys, ytot))';
  xpnames = cat(1, xpnames, xpnametot);
end
