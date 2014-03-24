function [s,H]=singularvalues_mechanism(km,x,y)
%  [s,H]=singularvalues_mechanism(km, x, y)
%
%  Returns the singular values of the linearized observation
%  function.
%
%  Input
%     km        ->   kinematic model struct. See build_model
%     x         ->   state vector
%     y         ->   Optional. observation vector. Useful for
%                    studying the effect of missing markers
%  Output
%     s         <-   vector with singular values
%     H         <-   Linearized observation function

% Kjartan Halvorsen
% 2004-03-22

% Prepare marker positions.

[initnames, p0, p0vec] = prepare_mdata(km.p0);

if nargin < 3
  y = ones(size(p0vec));
end
if nargin < 2
  x = zeros(length(km.gcnames),1);
end

[yy,H,Gs]=observe_mechanism_H(x,y,km.twists,p0);

H=H(1:length(y),1:length(x));

s=svd(H);
