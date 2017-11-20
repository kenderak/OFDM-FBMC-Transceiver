clear all;close all;clc;

%% OFDM/FBMC System Simulation - 2017/18/1
% Settings
simulationMethod = 'FBMC';   % OFDM or FBMC system simulation 
modulationMethod = '4QAM';   % BPSK,QPSK,4QAM,16QAM,64QAM
codingTechnique = 'None';    % None, ...
numOfSym = 1000;              % Number of symbols
sizeOfFFT = 128;              % Size of IFFT/FFT
numOfCarrier = 88;           % Number of data carriers 
overSampling = 2;            % Factor of oversampling (1,2,4 ...)
cpLength = 0;                % Cyclic prefix length for an OFDM symbol
K = 4;                       % Overlapping factor for FMBC modulation
CR = 7;                     % Clipping Ratio [dB]

% Set the parameter object
switch simulationMethod
    case 'OFDM'
        Param = paramOFDM(modulationMethod, numOfSym, sizeOfFFT, numOfCarrier, overSampling, cpLength);
    case 'FBMC'
        Param = paramFBMC(modulationMethod, numOfSym, sizeOfFFT, numOfCarrier, overSampling, K);
end

%% I. OFDM/FBMC Transmitter
% 1. Generate binary input data (symbols)
binaryData = randi([0 1], Param.M*numOfSym*numOfCarrier, 1);
% 2. Channel encoding 
encodedData = encoder(binaryData,codingTechnique);
% 3. Symbol mapper
mappedData = Param.mapper(encodedData);
% 4. Perform OFDM/FBMC modulation
switch simulationMethod
    case 'OFDM'
        Modulated = modulatorOFDM(mappedData,Param);
    case 'FBMC'
        Modulated = modulatorFBMC(mappedData,Param);
        %Modulated = modulatorFBMC_PPN(mappedData,Param);
        %Modulated = modulatorFBMC_PPN_NFFT(mappedData,Param);
        
        %Quantization
        %ModulatedQ = quant_single_modulatorFBMC(mappedData,Param);
end

%% Clipping
signal = Modulated.signalTx;
gamma = 10^(CR/20);
sigma = sqrt(Modulated.Es);
Amax = sigma*gamma;
phase = exp(1i*angle(signal));
signal(abs(signal) > Amax) = Amax*phase(abs(signal) > Amax);
%Modulated.signalTx = signal;

%% II. Settings of noise parameters for channel
% Ratio for BER curve
EbNo = 0:1:13; % Bit to Noise Ratio [dB]
EsNo = EbNo  + 10*log10(sizeOfFFT/(sizeOfFFT+Param.CP))+ 10*log10(Param.M*numOfCarrier/sizeOfFFT); % Symbol to Noise Ratio [dB]
SNR = EsNo; % Signal to Noise Ratio for the AWGN channel [dB]
theoryBer =  (1/2)*erfc(sqrt(10.^(EbNo/10)));
% Pre-allocate vectors for the errors
numOfErrors = zeros(1,length(SNR));
BER = zeros(1,length(SNR));

%% III. OFDM/FBMC Receiver
for i=1:length(SNR)
    % 0. Passing through the channel
    signalRx = awgn(Modulated.signalTx, SNR(i),'measured');
    % 1. Perform OFDM/FBMC demodulation
    switch simulationMethod
        case 'OFDM'
            Demodulated = demodulatorOFDM(signalRx, Param);
        case 'FBMC'
            Demodulated = demodulatorFBMC(signalRx, Param);
            %Demodulated = demodulatorFBMC_THETA(signalRx, Param);
    end
    % 2. Symbol demapping
    demappedData = Param.demapper(Demodulated.Symbols(:));
    % 3. Channel decoding
    decodedData = decoder(demappedData, codingTechnique);
    % Calculate the errors
    [numOfErrors(i),BER(i)]=biterr(binaryData,decodedData);
end
%% Performance Analysis
% Power Spectral Density
figure(1)
pwelch(Modulated.signalTx*sqrt(sizeOfFFT), hann(overSampling*sizeOfFFT),...
    [],overSampling*sizeOfFFT,overSampling,'centered'); hold on;

% pwelch(ModulatedQ.signalTx*single(sqrt(sizeOfFFT)), single(hann(overSampling*sizeOfFFT)),...
%     [],single(overSampling*sizeOfFFT),single(overSampling),'centered');

% Complementary Cumulative Distribution Function (CCDF) estimation
CCDF = ccdf(Modulated,Param,simulationMethod);

% Bit Error (BER) curves 
figure(3)
semilogy(EbNo,theoryBer,'bs-','LineWidth',1.5);
hold on;
semilogy(EbNo,BER,'ro--','linewidth',1.5);
grid on
axis([0 10 10^-6 1])
legend('BPSK theory ', [modulationMethod ,' simulation'])
xlabel('E_{b}/N_{o} [dB]')
ylabel('Bit Error Rate')
title(['Bit error probability curve for ', modulationMethod,  ' using ', simulationMethod])