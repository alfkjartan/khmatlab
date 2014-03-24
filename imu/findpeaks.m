function pks = findpeaks(x, cutoff, thr, peakwidth)

  [b2,a2] = butter(4, cutoff);
  xf = filtfilt(b2,a2,x);
  pks = peakdet(xf, 1e-3);

  abovethr = find(xf > thr);
  pks = intersect(pks, abovethr); 
  
  TR = pks(2:end) - pks(1:end-1);
  doubletsR = find(TR < peakwidth);  % Find repated pks. Keep the last
  pks(doubletsR) = [];
