function [jc, jax] = jc_estim_symmetry(leftmdata, rightmdata,... 
				       proximalmarkers, distalmarkers, refdata)
%  [jc, jax] = jc_estim_symmetry(mdata, proximalmarkers, ...
%                                distalmarkers, refdata)
%  Returns an estimate of the joint center and axis of rotation.
%  The marker data must contain sufficient movement of the distal limb.
%  Uses symmetry to estimate one jc for two joints, assumed to be 
%  positioned symmetric with respect to a plane defined by the 
%  proximal markers. Mostly used for estimating hjc. The proximal
%  markers should be {'asis_l', 'asis_r', 'spine'}, in that order.
%
% Input
%    mdata           ->   Cell array with marker data. 
%                         Each cell contains {attributes, data}
%    proximalmarkers ->   Cell array of marker names
%    distalmarkers   ->   Cell array with cell array of marker
%                         names {leftsegm, rightsegm}
%    refdata         ->   Optional. Marker data, {attributes, data},
%                         gives the position of the proximalmarkers in
%                         their reference position.
% Output
%    hjc            <-   the hip joint centers with respect to a
%                        coordinate system fixed to the proximal. If
%                        refdata is not provided, then the coordinate
%                        system coincides with the labsystem for
%                        the first frame of data. If refdata is
%                        provided, then the coordinate system
%                        coincides with the labsystem when the
%                        proximal markers are in the reference position

% Kjartan Halvorsen
% 2003-08-13

try
  
leftdistalrawdata = [];
leftproximaldata = [];
for i=1:length(leftmdata)
  lmnames = getvalue(leftmdata{i}{1}, 'MARKER_NAMES');
  proximaldata = extractmarkers(leftmdata{i}{2}, ...
			       lmnames, proximalmarkers);
  distalrawdata = extractmarkers(leftmdata{i}{2}, ...
				 lmnames, distalmarkers{1});
  leftdistalrawdata =  cat(1, leftdistalrawdata, distalrawdata);
  leftproximaldata =  cat(1, leftproximaldata, proximaldata);
end

rightdistalrawdata = [];
rightproximaldata = [];
for i=1:length(rightmdata)
  rmnames = getvalue(rightmdata{i}{1}, 'MARKER_NAMES');
  proximaldata = extractmarkers(rightmdata{i}{2}, ...
				rmnames, proximalmarkers);
  distalrawdata = extractmarkers(rightmdata{i}{2}, ...
				 rmnames, distalmarkers{2});
  try
  rightdistalrawdata =  cat(1, rightdistalrawdata, distalrawdata);
  rightproximaldata =  cat(1, rightproximaldata, proximaldata);
  catch
    keyboard
  end
end


% Load reference position data
proximalrefdata = extractmeanmarkers(refdata, ...
				 proximalmarkers);
proximalrefdata = (proximalrefdata(:))';

leftdistalref = extractmeanmarkers(refdata, ...
			       distalmarkers{1});
leftdistalref = (leftdistalref(:))';

rightdistalref = extractmeanmarkers(refdata, ...
			       distalmarkers{2});
rightdistalref = (rightdistalref(:))';

leftproximaldata = cat(1, proximalrefdata, leftproximaldata);
rightproximaldata = cat(1, proximalrefdata, rightproximaldata);
leftdistalrawdata = cat(1, leftdistalref, leftdistalrawdata);
rightdistalrawdata = cat(1, rightdistalref, rightdistalrawdata);

  
% Compute the coordinates of the distal markers in a proximal fixed
% coordinate system
[leftdistaldata,leftres] = getRelMotion(leftdistalrawdata, leftproximaldata);
[rightdistaldata,rightres] = getRelMotion(rightdistalrawdata, ...
					  rightproximaldata); 

residualthr = 10;
leftdistaldata = leftdistaldata(find(leftres<residualthr),:);
rightdistaldata = rightdistaldata(find(rightres<residualthr),:);

% Reflect the right side marker data
asis_l = proximalrefdata(1:3)';
asis_r = proximalrefdata(4:6)';
%spine = proximalrefdata(7:9)';
spine = mean([asis_l asis_r],2); % Use this when spine marker
                                  %position may be assymetric.
vertical = [0;0;1];

nvec = asis_l-asis_r;
nvec = nvec - (nvec'*vertical)*vertical;
nvec = nvec/norm(nvec);

sp = 2*nvec*nvec'*spine;
Refl = eye(3) - 2*nvec*nvec';

nrmrks = size(rightdistaldata,2)/3;
spp = kron(sp,ones(1,nrmrks));
for fr = 1:size(rightdistaldata,1)
  p = reshape(rightdistaldata(fr,:),3,nrmrks);
  if (sum(p==0) ~= 3*nrmrks)
    p = Refl*p + spp;
    rightdistaldata(fr,:) = reshape(p,1,3*nrmrks);
  end
end

distaldata = cat(1, leftdistaldata, rightdistaldata);

[nfr, nmrks] = size(distaldata);
nmrks = nmrks/3;

%keyboard

% Reshape and shift dimensions so that the marker data can be used
% in the biascomp_ls function
td = permute(reshape(distaldata', 3, nmrks, nfr), [1 3 2]);
[jc, jax] = biascomp_ls(td);
%hjc = hiniduma_ls(td);


% Reflect jc to get right side jc.
jc = cat(2, jc, Refl*jc + sp);

  
catch
  disp('ERROR OCCURRED IN JC_ESTIM_SYMMETRY')
  keyboard
end
