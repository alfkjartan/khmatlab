function drawwhiskerplot(draw_data, type, lineWidth, width, legnd)
% drawwhiskerplot(draw_data, type, lineWidth, width)
% 
% Input
%    draw_data   ->  n x 2 x m, or n x 3 x m,  where n is length of series, m is
%                    number of series
%     type:
%   0   - default. whiskers both ways
%   1   - upper whiskers only
%   2   - lower whiskers only

% Kjartan Halvorsen


symbols = {'ko', 'k+', 'ks', 'kv', 'kd'};

n = size(draw_data, 1);
m = size(draw_data, 3);

if nargin < 3
  width = 0.4 / m ;
  lineWidth = 1;
end

if nargin < 2
  type = 0;
end


if (size(draw_data, 2) == 2)
  unit = (1-1/(1+n))/(1+9/(width+3));
else
  range = max(draw_data(:,1,:)) - min(draw_data(:,1,:));
  unit = (1-1/(1+range))/(1+9/(width+3));
end

hold on;       

plhndls = zeros(1,m);
for j = 1:m
  offs = (j-ceil(m/2))*unit/3;
  for i = 1:n
        
    if (size(draw_data, 2) == 2)
      v = draw_data(i,:,j);
      xx = i;
    else
      v = draw_data(i,2:3,j);
      xx = draw_data(i,1,j);
    end
    
    vup = v(1)+v(2);
    vdn = v(1)-v(2);
    if type==1
      vd = vup;
    elseif type==2
      vup = vd;
    end
    
    % draw the mean
    plhndsl(j) = plot([xx+offs], v(1), symbols{j},  'LineWidth', lineWidth);
    
    if (type==0 | type == 1)
      % draw the max line
      plot([xx+offs-unit, xx+offs+unit], [vup, vup], 'LineWidth', lineWidth);
      % draw vertical line
      plot([xx+offs, xx+offs], [v(1), vup], 'LineWidth', lineWidth);
    end
    
    if (type==0 | type == 2)
      % draw min line
      plot([xx+offs-unit, xx+offs+unit], [vdn, vdn], 'LineWidth', lineWidth);
      plot([xx+offs, xx+offs], [v(1), vdn], 'LineWidth', lineWidth);
    end
    
  end
end

if (nargin > 4)
  legend(plhndsl, legnd)
end
