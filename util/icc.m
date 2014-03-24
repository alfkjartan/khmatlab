function ic = icc(d1,d2)
  %  ic = icc(d1,d2)
  % Computes the intra-class correlation coefficient.
    
  % Kjartan Halvorsen
  % 2004-05-14
   
  d1=d1(:);
  d2=d2(:);

  totmean = mean(cat(1,d1,d2));
  rowmeans = mean(cat(2,d1,d2),2);
  
  SSr = 2* sum((rowmeans - totmean).^2);
  
  SSw = sum((cat(1,d1-rowmeans, d2-rowmeans)).^2);
  
  n=length(d1);
  N = 2*n;
  MSr = SSr/(n-1);
  MSw = SSw/(N-n);
  
%  ic = (MSr - MSw)/(MSr);
  ic = (MSr - MSw)/(MSr);
  
  if 0 % Old stuff
    stdd = std(cat(2,d1,d2),0,2);
  
    withinsubjvar = mean(stdd.^2);
  
    betweensubjvar = var(mean(cat(2,d1,d2),2));
  
    ic = betweensubjvar / (betweensubjvar + withinsubjvar);
  
  end
  