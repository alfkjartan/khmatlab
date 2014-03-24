function [data_n, h0] = split_and_present(data, splitcolumn, dt, lgnds)
%  [data_n, figh] = split_and_present(data, splitcolumn, dt, lgnds)
%
% Splits a data set with cyclical data into separate cycles,
% normalizes to 100 data points per cycle, and presents the data as
% mean +/- stdv.
%
% Input
%    data        ->  data matrix (nfr x nchs)
%    splitcolumn ->  column of data for which to find start and end
%                    of cycles.
%    dt          ->  sampling interval. Optional
%    lgnds       ->  legends. Cell array with strings. Optional 
% Output
%    data_n      <-  normalizes data matrix (100 x ncycles x nchs)
%    figh        <-  handle to the figure
%


% Kjartan Halvorsen
% 2004-10-25

if (nargin < 3)
  dt = 1;
end

samplefreq = 1/dt;

nchs = size(data,2);

% ----------------------------------
% Identify start and end of cycles
% ----------------------------------

[cyclfr] = detectcycles(data(:, splitcolumn), ...
			samplefreq, samplefreq, 'displacement');

[cyclefr,figh]=detectstepsplot(cyclfr, ...
			       data(:, splitcolumn));
close(figh);
  
cyclefr(1)=max(cyclefr(1),1);

% Remove repetitions in cyclefr
cyclefr=unique(cyclefr);
ncycles=length(cyclefr)-1;

% ----------------------------------
% Determine duration of each cycle
% ----------------------------------

nsmp=zeros(ncycles,1);
for i = 1:ncycles
  nsmp(i)=cyclefr(i+1)-cyclefr(i);
end
  
cycledur=nsmp*(1/samplefreq);


% -----------------------------------------
% Divide the data according to the cycles, 
% then normalize to 100 data points
% -----------------------------------------
data_n = zeros(100,ncycles, nchs);
for j=1:nchs
  for c=1:ncycles
    data_n(:,c,j) = interp1((1:nsmp(c))'*100/nsmp(c), ...
			    data(cyclefr(c):(cyclefr(c+1)-1),j), ...
			    (1:100)');
  end
end


% -----------------------------------------
% Present results
% -----------------------------------------

h0 = figure('Color',[1 1 1], ...
	    'Name', ['cyclic data'], ...
	    'NumberTitle','off', ...
	    'PaperPosition',[18 180 573 432], ...
	    'PaperType', 'a4',...
	    'PaperUnits','points', ...
	    'Position',[20 20 600 610], ...
	    'ToolBar','none');


plotcolororder = {'b','m','g', 'r', 'b--', 'm--', 'g--', 'r--'};



angle_mean = mean(data_n,2);  
angle_std = std(data_n,0,2);  
lh=[];
for j=1:nchs
  pclr= plotcolororder{mod(j,length(plotcolororder)/2+1)};
  lhj = plot((-10:110), ...
	    cat(1, angle_mean(90:100,1,j), ...
		angle_mean(:,1,j), ...
		angle_mean(1:10,1,j)), ...
	    pclr,...
	    (-10:110), ...
	    cat(1, angle_mean(90:100,1,j) + angle_std(90:100,1,j), ...
		angle_mean(:,1,j) + angle_std(:,1,j), ...
		angle_mean(1:10,1,j) + angle_std(1:10,1,j)), ...
	    [pclr,':'],...
	    (-10:110), ...
	    cat(1, angle_mean(90:100,1,j) - angle_std(90:100,1,j), ...
		angle_mean(:,1,j) - angle_std(:,1,j), ...
		angle_mean(1:10,1,j) - angle_std(1:10,1,j)), ...
	    [pclr,':']);
  hold on
  lh = cat(1, lh, lhj);
end

ylim = get(gca,'YLim');
xlim = [-10 110];
set(gca,'XLim', xlim);

plot([0 0], ylim, 'k:');
plot([100 100], ylim, 'k:');
plot(xlim, [0 0], 'k:');
  
title('Normalized data', 'FontSize', 12)
xlabel('% of cycle', 'FontSize', 12);


if (nargin > 3)
  legend(lh(1:3:end), lgnds{:})
end

  