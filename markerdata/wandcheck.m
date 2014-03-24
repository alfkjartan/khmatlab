%function wandcheck(filename, wandlength)
  % Computes the distance between the two wandmarkers, and compare
  % to the known length of the wand. 
  % The discrepancy gives a measure of the accuracy of the
  % setup. The result is presented at orthogonal cross sections of the of
  % measurement volume through its center.
  
  % Kjartan Halvorsen
  % 2003-10-27
    
  fid = fopen(filename);
  
  if (fid == -1)
    error(['Could not open file ', filename])
  end
  
  [attr, md] = read3dtsv(fid);
  
  fclose(fid);
  
  nfr = size(md,1);
  
  res = zeros(nfr,4);

  removeind = [];
  for fr = 1:nfr
    m1 = md(fr,1:3);
    m2 = md(fr, 4:6);
    
    if (~m1 | ~m2)
      removeind = cat(1, removeind, fr);
    end
    
    res(fr,1:3) = 0.5*m1 + 0.5*m2;
    
    res(fr,4) = norm(m1-m2) - wandlength;
  end

  res(removeind,:) = [];
  cdepth = 128;
  cmap = jet(cdepth);
  
  abserr = abs(res(:,4)); 
  minerr = min(abserr);
  erange = max(abserr) - minerr;

  figure
  for fr = 1:nfr
    cind = fix((abserr(fr) - minerr)/erange*(cdepth-1)) + 1;
    plot3(res(fr,1), res(fr,2), res(fr,3), 'o', ...
	  'MarkerEdgeColor', cmap(cind,:));
    hold on
  end
  
