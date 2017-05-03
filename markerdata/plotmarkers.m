function [h, sqerr]=plotmarkers(varargin)
%  [h, sqerr]=plotmarkers(md1,md2) or
%  [h, sqerr]=plotmarkers(md1,mnames1,md2,mnames2,figs,lgnd)
% 
% Plots the error in the marker data, each marker in a separate figure 
% window, and the coordinates in separate graphs.
%
% Input: 
%     md1,md2     ->   matrices (nfr x nm*3) with marker data.
%     mnames      ->   (optional) cell array (size nm) of strings.
%     figs        ->   (optional) vector of figure handles, or nx2 cell
%                      array where the first column contains names
%                      of the markers to plot, and the second
%                      column corresponding figure handles 
%     lgnd        ->   (optional) strings for the legend
% Output
%     h           <-   handles to the figures
%     sqerr       <-   mean (over the set of markers) of squared
%                      errors for each frame.
  

% Kjartan Halvorsen
% 2001-07-24

% Revisions
% 2003-11-18   Using varargin, and possibly two arguments only
% 2003-12-08   Output squared errors.
% 2009-06-29   Added optional cell array with names and figure
%              handles as input argument number five.


%%try
  
md1 = varargin{1};
[nfrs,cols]=size(md1);
nms = cols/3;

if (nargin<6)
  legnd={'series 1','series 2'};
else
  legnd = varargin{6};
end
if (nargin < 5)
  figs = zeros(nms,1);
  for f=1:nms
    figs(f) = figure;
  end
else
  if iscell(varargin{5})
    m2plot = varargin{5};
    mnames = varargin{2};
    md1 = extractmarkers(md1,mnames,m2plot(:,1));
    varargin{3} = extractmarkers(varargin{3},varargin{4},m2plot(:,1));
    varargin{2} = m2plot(:,1);
    varargin{4} = m2plot(:,1);
    figs = cat(1, m2plot{:,2});
  end
end
if (nargin == 2)
   md2 = varargin{2};
   mnames1 = cell(nms,1);
   for m=1:nms
     mnames1{m} = sprintf('marker %d', m);
   end
   mnames2 = mnames1;
else
  mnames1 = varargin{2};
  md2 = varargin{3};
  mnames2 = varargin{4};
end
  
[nfrs2,cols2]=size(md2);

T1=(1:nfrs);
T2=(1:nfrs2);

nfr = min(nfrs, nfrs2);
sqerr=zeros(nfr, nms);
for i=1:min(length(figs),nms)
   marker=mnames1{i};
   [slask,mind]=intersect(mnames2, {marker});

%   if 0
   for j=1:nfr
     if (~hasmissing(md1(j, (i-1)*3+1:i*3)) & ...
	 ~hasmissing(md2(j, (mind-1)*3+1:mind*3)))
       sqerr(j,i) = sqrt( ...
       (md1(j,(i-1)*3+1) - md2(j,(mind-1)*3+1)).^2 + ...
       (md1(j,(i-1)*3+2) - md2(j,(mind-1)*3+2)).^2 + ...
       (md1(j,(i-1)*3+3) - md2(j,(mind-1)*3+3)).^2 );
     end
   end
%   end

   %keyboard
   figure(figs(i))
   subplot(311)
   plot(T1,md1(:,(i-1)*3+1),'r',T2,md2(:,(mind-1)*3+1),'b')
   legend(legnd{:}, 'Location','best')
   subplot(312)
   plot(T1,md1(:,(i-1)*3+2),'r',T2,md2(:,(mind-1)*3+2),'b')
   legend(legnd{:}, 'Location', 'best')
   subplot(313)
   plot(T1,md1(:,(i-1)*3+3),'r',T2,md2(:,(mind-1)*3+3),'b')
   legend(legnd{:}, 'Location', 'best')
   
   ylabel('z')
   subplot(312)
   ylabel('y')
   subplot(311)
   ylabel('x')
   title(marker, 'Interpreter', 'none')
end

h=figs;

sqerrtmp=zeros(nfr,1);
for i=1:nfr
    msqi = mean(sqerr(i,find(sqerr(i,:) ~= 0)));
    if ~isempty(msqi)
      sqerrtmp(i) = msqi;
    end
end

sqerr=sqerrtmp;

%%catch 
%%  keyboard
%%end
