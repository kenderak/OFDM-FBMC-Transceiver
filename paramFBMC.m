function [ Param ] = paramFBMC( mod, numOfSym, sizeOfFFT, numOfCarrier, overSampling, K )

Param.N     = sizeOfFFT;
Param.D     = numOfCarrier;
Param.K     = K;
Param.S     = Param.N * Param.K;
Param.OV    = overSampling;
Param.Offset = Param.N/2;
Param.CP    = 0;

if (numOfCarrier == sizeOfFFT)
Param.CarrierIndexes    = [1:Param.N]';
else
Param.CarrierIndexes    = [2:Param.D/2+1,Param.N-Param.D/2+1:Param.N]';
end

switch mod
    case 'BPSK'
        Param.M = 1;
    case '4QAM'
        Param.M = 2;
    case '16QAM'
        Param.M = 4;
    case '64QAM'
        Param.M = 6;
    otherwise
        Param.M = 2;
end

Param.NrOfSymbols = numOfSym;
Param.NrOfBits = numOfCarrier*numOfSym*Param.M;

% Configurate symbol mapper/demapper
Param.mapper = comm.RectangularQAMModulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');
Param.demapper = comm.RectangularQAMDemodulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');

% Prototype filter coefficients with correct phase
switch Param.K
    case 4
    Param.H_coeffs   = ([1, 0.97196, sqrt(2)/2, .235147, 0 ,0.235147, sqrt(2)/2, 0.97196].*(-1).^([0:7]))';
    case 2
    Param.H_coeffs = ([1, sqrt(2)/2, 0 , sqrt(2)/2].*(-1).^([0:3]))';
end
% Build the filterbank
Param.H_coeffsFB = circshift(Param.H_coeffs, Param.K);

    Param.FB_odd       = kron(real(ones(1,Param.N/2)),Param.H_coeffsFB);
    Param.FB_odd       = Param.FB_odd(:);
    
    Param.FB_even      = kron(real(ones(1,Param.N/2)),Param.H_coeffsFB);
    Param.FB_even      = circshift(Param.FB_even(:),-4);

    Param.EsH = mean(abs(Param.H_coeffs).^2);

end

