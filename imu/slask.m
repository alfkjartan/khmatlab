figure(14)
clf
subplot(211)
%%plot(TA,rmse_repA, 'm', 'linewidth', 6)
plot(TA,rmseA, 'm', 'linewidth', 6)
hold on
%%plot(T,rmse_repB, 'r', 'linewidth', 6)
plot(T,rmseB, 'r', 'linewidth', 6)
%%plot(time1,rmse_repref, 'b', 'linewidth', 6)
plot(time1,rmse_ref, 'b', 'linewidth', 6)

legend('Algorithm 2G', 'Algorithm 2GA', 'LP filter', 'location', 'northeast')


ylabel('ang vel [rad/s]', 'fontsize', 16);
title('RMS error')

subplot(212)
plot(time1, wb1true', 'linewidth', 4);

ylabel('ang vel [rad/s]', 'fontsize', 16);
xlabel('time [s]', 'fontsize', 16);
title('True angular velocity')

set(findobj(gcf, '-property', 'fontsize'), 'fontsize', fntsze)
print(fullfile(respth, "RMSerror.eps"), "-depsc");

