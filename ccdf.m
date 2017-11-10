function [CCDF] = ccdf(Modulated,Param, sim)
% Calculate and plot CCDF curve
signalTx = Modulated.signalTx;
%% Parameters
S = Param.S;
N = Param.N;
OV = Param.OV;
CP = Param.CP;
Offset = Param.Offset;      

% Peak-to-Average Ratio (PAPR[dB]) estimation

switch sim
    case 'OFDM'
        
        numOfSym = Param.numOfSym;
        PAPR_dB = zeros(numOfSym,1);
        
        for i = 1:numOfSym
            index_start = 1 + (i-1)*(N+CP)*OV;
            index_end = i * (N+CP)*OV;
            
            PAPR_dB(i,1) = 10*log10(max(abs(signalTx(index_start:index_end))).^2 ...
                / mean(abs(signalTx(index_start:index_end)).^2));
        end
        
    case 'FBMC'
        
        numOfSym = Param.NrOfSymbols;
        PAPR_dB = zeros(numOfSym-19,1);
        
        for ii=10:numOfSym-10
            
            index_start = 1 + (ii-1)*N*OV;
            index_end   = ii * N*OV;
             
            PAPR_dB(ii,1) = 10*log10(max(abs(signalTx(index_start:index_end))).^2 ...
                / mean(abs(signalTx(index_start:index_end)).^2));
        end
end

% Plot CCDF curves
[ccdf_y,ccdf_x] = ecdf(PAPR_dB);
ccdf_y = 1 - ccdf_y;
ccdf_orig = 1-(1-exp(-10.^(ccdf_x/10))).^N;

CCDF.PAPRdB = PAPR_dB;
CCDF.ccdf_y = ccdf_y;
CCDF.ccdf_x = ccdf_x;
CCDF.ccdf_orig = ccdf_orig;

figure(2);
semilogy(ccdf_x,ccdf_y,'-r','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3),grid on;
hold on;
semilogy(ccdf_x,ccdf_orig,'-b','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3),grid on;
title('Complementary Cumulative Distribution Function of the PAPR')
xlabel('PAPR_{0} [dB]');
ylabel('CCDF(PAPR) = P(PAPR > PAPR_{0})');
legend([sim, ', N =' num2str(N)] , ['Original, N =' num2str(N)]);
xlim([0 14]); ylim ([10^-4 1])
