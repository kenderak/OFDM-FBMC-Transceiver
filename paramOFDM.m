function Param = paramOFDM(mod, numOfSym, N, D, OV, CP)

Param.N = N;
Param.D = D;
Param.OV = OV;
Param.CP = CP;

if (D == N)
Param.CarrierIndexes    = [1:Param.N]';
else
Param.CarrierIndexes    = [2:Param.D/2+1,Param.N-Param.D/2+1:Param.N]';
end

switch mod
    case 'BPSK'
        Param.M = 1;
        Param.mapper = comm.BPSKModulator;
        Param.demapper = comm.BPSKDemodulator;
    case 'QPSK'
        Param.M = 2;
        Param.mapper = comm.QPSKModulator('BitInput', true);
        Param.demapper = comm.QPSKDemodulator('BitOutput', true);
    case '4QAM'
        Param.M = 2;
    case '16QAM'
        Param.M = 4;
    case '64QAM'
        Param.M = 6;
end

switch mod
    case '4QAM'
        Param.mapper = comm.RectangularQAMModulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');
        Param.demapper = comm.RectangularQAMDemodulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
    case '16QAM'
        Param.mapper = comm.RectangularQAMModulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');
        Param.demapper = comm.RectangularQAMDemodulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
    case '64QAM'
        Param.mapper = comm.RectangularQAMModulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');
        Param.demapper = comm.RectangularQAMDemodulator(...
    'ModulationOrder', 2^Param.M, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
end

Param.numOfSym = numOfSym;
Param.numOfBits = D*Param.M*numOfSym;

end