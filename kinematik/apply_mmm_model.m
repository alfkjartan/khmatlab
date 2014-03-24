function [com, am] = apply_mmm_model(md, markers_com, weights_com, ...
				     markers_am, weights_am)
% [com, am] = apply_mmm_model(md, markers_com, weights_com, ...
%				     markers_am, weights_am)
%
% Main function for applying the "Momentum from Minimal set of
% Markers" method. 
%


% Kjartan Halvorsen
% 2008-01-09


ncom_m = size(markers_com,1);
found_m = 0;
col = 1;
while ( (found_m < ncom_m) & (col <= size(markers_com,2)) )
  [mm, slask, found_m] = extractmarkers(md, markers_com(:,col));
  col = col + 1;
  disp([' Markers found: ' , int2str(found_m)])

end


for i=1:length(weights_com)
  mm(:,(i-1)*3+1:i*3) =   mm(:,(i-1)*3+1:i*3) * weights_com(i);
end

com = cat(2, ...
	  sum(mm(:,1:3:end), 2), ...
	  sum(mm(:,2:3:end), 2), ...
	  sum(mm(:,3:3:end), 2));

if nargin > 3
  % Compute also angular momentum
  
  fs = str2double(getvalue(md{1}, 'FREQUENCY'));
  
  ncom_m = size(markers_am,1);
  found_m = 0;
  col = 1;
  while ( (found_m < ncom_m) & (col <= size(markers_am,2)) )
    [mm,slask,found_m] = extractmarkers(md, markers_am(:,col));
    col = col + 1;
  disp([' Markers found: ' , int2str(found_m)])
  end
  
  vel = centraldiff(mm,fs);

  [nr,nc] = size(mm);
  
  mmxveltmp = cross(cat(2, ...
			reshape(mm(:,1:3:end), nr*nc/3, 1),...
			reshape(mm(:,2:3:end), nr*nc/3, 1),...
			reshape(mm(:,3:3:end), nr*nc/3, 1)), ...
		    cat(2, ...
			reshape(vel(:,1:3:end), nr*nc/3, 1),...
			reshape(vel(:,2:3:end), nr*nc/3, 1),...
			reshape(vel(:,3:3:end), nr*nc/3, 1)));
  mmxvel = zeros(size(mm));
  mmxvel(:,1:3:end) = reshape(mmxveltmp(:,1), nr, nc/3);
  mmxvel(:,2:3:end) = reshape(mmxveltmp(:,2), nr, nc/3);
  mmxvel(:,3:3:end) = reshape(mmxveltmp(:,3), nr, nc/3);
  
for i=1:length(weights_am)
  mmxvel(:,(i-1)*3+1:i*3) =   mmxvel(:,(i-1)*3+1:i*3) * weights_am(i);
end

% Angular momentum wrt origo  
am_wrt_origo = cat(2, ...
	  sum(mmxvel(:,1:3:end), 2), ...
	  sum(mmxvel(:,2:3:end), 2), ...
	  sum(mmxvel(:,3:3:end), 2));
  
% Subtract angular momentum of the com
comvel = centraldiff(com, fs);
com_am = cross(com,comvel);
am = am_wrt_origo - com_am;
%am = am_wrt_origo;

%keyboard
  
end

