function [jc, jax] = jc_estim(mdata, proximalmarkers, distalmarkers, refdata)
%  [jc, jax] = jc_estim(mdata, proximalmarkers, distalmarkers, refdata)
%  Returns an estimate of the joint center and axis of rotation.
%  The marker data must contain sufficient movement of the distal limb.

% Input
%    mdata           ->   marker data. Cell array: {attributes, data}
%    proximalmarkers ->   Cell array of marker names
%    distalmarkers   ->   Cell array of marker names
%    refdata         ->   Optional. Marker data, {attributes, data},
%                         gives the position of the proximalmarkers in
%                         their reference position.
% Output
%    hjc            <-   the hip joint center with respect to a
%                        coordinate system fixed to the proximal. If
%                        refdata is not provided, then the coordinate
%                        system coincides with the labsystem for
%                        the first frame of data. If refdata is
%                        provided, then the coordinate system
%                        coincides with the labsystem when the
%                        proximal markers are in the reference position

% Kjartan Halvorsen
% 2003-07-29
  
distaldata = extractmarkers(mdata, distalmarkers);

if ~isempty(proximalmarkers)
  proximaldata = extractmarkers(mdata, proximalmarkers);

  % Compute the coordinates of the distal markers in a proximal fixed
  % coordinate system
  distaldata = getRelMotion(distaldata, proximaldata);
end

[nfr, nmrks] = size(distaldata);
nmrks = nmrks/3;

% Reshape and shift dimensions so that the marker data can be used
% in the biascomp_ls function
td = permute(reshape(distaldata', 3, nmrks, nfr), [1 3 2]);
[jc, jax] = biascomp_ls(td)
%hjc = hiniduma_ls(td);

if (nargin > 3)
  proximalrefdata = extractmarkers(refdata, proximalmarkers);
%  proximalrefdata = clipmdata(proximalrefdata);
  proximalrefdata = mean(proximalrefdata);
  
  T = soder(cat(1,proximaldata(1,:), proximalrefdata));
  
  jc = T * cat(1, jc, 1);
  jc = jc(1:3);
  jax = T * cat(1, jax, 0);
  jax = jax(1:3);
end


  