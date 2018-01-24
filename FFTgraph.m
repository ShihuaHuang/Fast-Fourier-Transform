plot(GOOG_Sept_20_2016,'g')
xlabel('Seconds');
ylabel('GOOG Price 20 September 2016');
title('The original, and filtered-price-paths for GOOG on 9/20/2016');
hold
plot(GOOG_5terms,'r')
plot(GOOG_10terms,'k')
plot(GOOG_15terms,'b')
legend('GOOG','GOOG (5 terms)','GOOG (10 terms)','GOOG (15 terms)')

plot(MSFT_Sept_20_2016,'g')
xlabel('Seconds');
ylabel('MSFT Price 20 September 2016');
title('The original, and filtered-price-paths for MSFT on 9/20/2016');
hold
plot(MSFT_5terms,'r')
plot(MSFT_10terms,'k')
plot(MSFT_15terms,'b')
legend('MSFT','MSFT (5 terms)','MSFT (10 terms)','MSFT (15 terms)')