function [] = ccdf(signalTx,sim,N)
% Calculate and plot CCDF curve

% Peak-to-Average Ratio (PAPR[dB]) estimation
numOfSym = size(signalTx,1);
PAPR_dB = zeros(1,numOfSym);


for ii=1:length(PAPR_dB)
    PAPR_dB(ii) = 10*log10(max(abs(signalTx(ii,:))).^2 / mean(abs(signalTx(ii,:)).^2));
end

% Plot CCDF curves
[ccdf_y,ccdf_x] = ecdf(PAPR_dB);
ccdf_y = 1 - ccdf_y;
ccdf_orig = 1-(1-exp(-10.^(ccdf_x/10))).^N;

figure(2);
semilogy(ccdf_x,ccdf_y,'-r','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3),grid on;
hold on;
semilogy(ccdf_x,ccdf_orig,'-b','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3),grid on;
title('Complementary Cumulative Distribution Function of the PAPR')
xlabel('PAPR_{0} [dB]');
ylabel('CCDF(PAPR) = P(PAPR > PAPR_{0})');
legend([sim, ', N =' num2str(N)] , ['Original, N =' num2str(N)]);
xlim([0 14]); ylim ([10^-4 1])
