function [ Modulated ] = quant_single_modulatorFBMC( ModulationSymbols, Param )

%% Parameters
S = single(Param.S);
N = single(Param.N);
OV = single(Param.OV);
Offset = single(Param.Offset);      
D = single(Param.D);
K = single(Param.K);
FB_odd = single(Param.FB_odd);
FB_even = single(Param.FB_even);
CarrierIndexes = single(Param.CarrierIndexes);

Modulated.NrOfSymbols = single(ceil(length(ModulationSymbols)/D));
Modulated.NrOfExtModSymbs = single(mod(D - mod(length(ModulationSymbols),D),D));
ModulationSymbols = single(ModulationSymbols);
%%
SymbolsF(CarrierIndexes,:) = single(reshape([ModulationSymbols;zeros(Modulated.NrOfExtModSymbs,1)],D,Modulated.NrOfSymbols));
Modulated.SymbolsF = SymbolsF;
odd_indexes  = CarrierIndexes(mod(CarrierIndexes,2)==0);
OddIndexesSpread= (odd_indexes-1)*K+1;
even_indexes = CarrierIndexes(mod(CarrierIndexes,2)==1);
EvenIndexesSpread=(even_indexes-1)*K+1;

% Normal symbols
SymbolsFOv(S,Modulated.NrOfSymbols) = single(eps*1i);
SymbolsFOv(OddIndexesSpread,:)= real(SymbolsF(odd_indexes,:));
SymbolsFOv(EvenIndexesSpread,:)= SymbolsFOv(EvenIndexesSpread,:)+ 1i*real(SymbolsF(even_indexes,:));
Modulated.SymbolsFOv = SymbolsFOv;

SymbolsFSpread = single(SymbolsFOv*0);
for i=-K:K
    SymbolsFSpread = SymbolsFSpread + circshift(SymbolsFOv,i);
end
Modulated.SymbolsFSpread = SymbolsFSpread;

% OFFSETSYMBOLS
SymbolsFOvOff(S,Modulated.NrOfSymbols)=single(eps*1i);
SymbolsFOvOff(OddIndexesSpread,:)= 1i*imag(SymbolsF(odd_indexes,:));
SymbolsFOvOff(EvenIndexesSpread,:)= SymbolsFOvOff(EvenIndexesSpread,:)+ imag(SymbolsF(even_indexes,:));
Modulated.SymbolsFOvOff = SymbolsFOvOff;

SymbolsFSpreadOff = single(SymbolsFOvOff*0);
for i=-K:K
    SymbolsFSpreadOff = SymbolsFSpreadOff + circshift(SymbolsFOvOff,i);
end
Modulated.SymbolsFSpreadOff = SymbolsFSpreadOff;

% Filtering in frequency domain
SymbolsFSpreadFilt(N,Modulated.NrOfSymbols) = eps*i;
SymbolsFSpreadOffFilt(N,Modulated.NrOfSymbols) = eps*i;

for i=1:Modulated.NrOfSymbols
    SymbolsFSpreadFilt(:,i) = real(Modulated.SymbolsFSpread(:,i)).*FB_odd +...
                                        1i*imag(Modulated.SymbolsFSpread(:,i)).*FB_even;
    SymbolsFSpreadOffFilt(:,i)= 1i*imag(Modulated.SymbolsFSpreadOff(:,i)).*FB_odd+...
                                           real(Modulated.SymbolsFSpreadOff(:,i)).*FB_even;
    
end
Modulated.SymbolsFSpreadFilt = SymbolsFSpreadFilt;
Modulated.SymbolsFSpreadOffFilt = SymbolsFSpreadOffFilt;

%% Time domain
Scale = single(1);
% Oversampling
SymbolsTOv = zeros(S * OV,Modulated.NrOfSymbols);
SymbolsTOv (1:S/2,:)= SymbolsFSpreadFilt(1:S/2,:);
SymbolsTOv(end-S/2+1:end,:)= SymbolsFSpreadFilt(S/2+1:end,:);
SymbolsTOv = ifft(SymbolsTOv) * Scale * OV;
Modulated.SymbolsTOv = SymbolsTOv;

SymbolsTOffOv = zeros(S*OV,Modulated.NrOfSymbols);
SymbolsTOffOv (1:S/2,:)= SymbolsFSpreadOffFilt(1:S/2,:);
SymbolsTOffOv (end-S/2+1:end,:)= SymbolsFSpreadOffFilt(S/2+1:end,:);
SymbolsTOffOv = ifft(SymbolsTOffOv) * Scale*OV;
Modulated.SymbolsTOffOv = SymbolsTOffOv;

Modulated.Scale = Scale;

% signal memory allocation;
signalTx = single(zeros((N/2 + N * Modulated.NrOfSymbols + (K-1)*N)*OV,1));
for i=1:Modulated.NrOfSymbols
    
    index_start = 1 + (i-1)*N*OV;
    index_end   =  (i-1) * N*OV + S*OV;
        
    % Normal symbols
    signalTx(index_start : index_end) = signalTx(index_start : index_end)+ SymbolsTOv(:,i);
    % Offseted symbols - complex values
    signalTx(index_start + Offset * OV : index_end + Offset * OV) = signalTx(index_start + Offset*OV:index_end + Offset*OV) + SymbolsTOffOv(:,i);

    
end
Modulated.signalTx = signalTx;
Modulated.Es = mean(abs(Modulated.signalTx).^2);
end
