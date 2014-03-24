  figure(7)
  clf
  subplot(211)
  plot(TB, xB(13:15,:)')
  subplot(212)
  plot(TB, q12B')
  hold on
  plot(TA, q12A', 'linewidth', 3)
