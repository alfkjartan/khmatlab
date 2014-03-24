function [p,r,delt] = linefit(varargin)
%  [a,r] = linefit(x, y [, xlabel, ylabel, title, legend])
%
% Plots x vs y and fits a first order polynomial (straight
% line). Parameters returned in a. 
% Calculates the r value (correlation coefficient).
% x and y may be cell arrays with different sets of data

% Kjartan Halvorsen
% 2004-02-19

x=varargin{1};
y=varargin{2};

if isnumeric(x);
  x = {x};
end
if isnumeric(y);
  y = {y};
end

if nargin > 2
  xlbl = varargin{3};
else
  xlbl = 'x';
end

if nargin > 3
  ylbl = varargin{4};
else
  ylbl = 'y';
end

if nargin > 4
  ttle = varargin{5};
else
  ttle = '';
end

if nargin > 5
  legnd = varargin{6};
else
  legnd = [];
end

if (length(x) ~= length(y))
  error('Data sets must have the same length')
end

markerstyles = {'o','+','s','x','d','v','*','h','.','<','>', ...
		'^'};

nstyles = length(markerstyles);

for i=1:length(x)
  x{i} = x{i}(:);
  y{i} = y{i}(:);
%  xpl = sort(unique(x{i}));
%  ypl = sort(unique(y{i}));
  xpl = x{i};
  ypl = y{i};
  
  if (length(xpl)>20)
    xpl=xpl(1:2:end);
    ypl=ypl(1:2:end);
  end
  
    
  mi = mod(i, nstyles);
  
  plot(xpl, ypl, markerstyles{mi})
  hold on
  
end


ylabel(ylbl);
xlabel(xlbl);
title(ttle,'Interpreter','none');

if ~isempty(legnd)
  legh=legend(legnd,0);
  chh = get(legh,'Children');
  for j=1:length(chh)
    try
      set(chh(j),'Interpreter','none');
    end
  end
end

%xlim = get(gca,'XLim');
%ylim = get(gca,'YLim');

%set(gca, 'XLim', [0 xlim(2)]);
%set(gca, 'YLim', [0 ylim(2)]);


xall = cat(1, x{:});
yall = cat(1, y{:});

[p, ps] = polyfit(xall, yall, 1);

[xx,inds] = sort(xall);
yall_sorted = yall(inds);

[yy, delt] = polyval( p, xx, ps);

hold on

plot(xx, yy, 'r','Linewidth',1.5)
plot(xx, yy+2*delt, 'r:','LineWidth',1.5)
plot(xx, yy-2*delt, 'r:','LineWidth',1.5)


text(0.53*max(xx)+0.5*min(xx), ...
     0.2*max(yall) + 0.8*min(yall),...
     ['y = ', num2str(p(1)), 'x + ', num2str(p(2))]);

% Compute r-value (correlation coefficient)
r=corrcoef(xall,yall);
  
% Do F-test
ssm = sum( (yy - mean(yall)).^2);
sse = sum( (yall_sorted - yy).^2);
sst = sum( (yall_sorted - mean(yall)).^2);

r2 = ssm/sst;
msm = ssm;
mse = sse / ( length(yall) - 2);

F = msm/mse;


text(0.53*max(xx)+0.5*min(xx), ...
     0.1*max(yall) + 0.9*min(yall),...
     ['r = ', num2str(r(1,2))])
text(0.53*max(xx)+0.5*min(xx), ...
     0.0*max(yall) + 1*min(yall),...
     ['r^2 = ', num2str(r2), '  F= ', num2str(F)])

nmin2 = length(yall)-2


  