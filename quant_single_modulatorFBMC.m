function Modulated = quant_single_modulatorFBMC(ModulationSymbols, Param)

%% Parameters
S = Param.S;
N = Param.N;
OV = single(Param.OV);
Offset = Param.Offset;      
D = Param.D;
K = Param.K;
CarrierIndexes = Param.CarrierIndexes;
FB_odd = single(Param.FB_odd);
FB_even = single(Param.FB_even);

NrOfSymbols = ceil(length(ModulationSymbols)/D);
NrOfExtModSymbs = mod(D - mod(length(ModulationSymbols),D),D);

%%
SymbolsF(CarrierIndexes,:) = single(reshape([ModulationSymbols;zeros(NrOfExtModSymbs,1)],D,NrOfSymbols));

odd_indexes  = CarrierIndexes(mod(CarrierIndexes,2)==0);
OddIndexesSpread= (odd_indexes-1)*K+1;
even_indexes = CarrierIndexes(mod(CarrierIndexes,2)==1);
EvenIndexesSpread=(even_indexes-1)*K+1;

% Normal symbols
SymbolsFOv(S,NrOfSymbols) = single(eps*1i);
SymbolsFOv(OddIndexesSpread,:)= real(SymbolsF(odd_indexes,:));
SymbolsFOv(EvenIndexesSpread,:)= SymbolsFOv(EvenIndexesSpread,:)+ 1i*real(SymbolsF(even_indexes,:));

SymbolsFSpread = single(SymbolsFOv*0);
for i=-K:K
    SymbolsFSpread = SymbolsFSpread + circshift(SymbolsFOv,i);
end

% OFFSETSYMBOLS
SymbolsFOvOff(S,NrOfSymbols)=single(eps*1i);
SymbolsFOvOff(OddIndexesSpread,:)= 1i*imag(SymbolsF(odd_indexes,:));
SymbolsFOvOff(EvenIndexesSpread,:)= SymbolsFOvOff(EvenIndexesSpread,:)+ imag(SymbolsF(even_indexes,:));

SymbolsFSpreadOff = single(SymbolsFOvOff*0);
for i=-K:K
    SymbolsFSpreadOff = SymbolsFSpreadOff + circshift(SymbolsFOvOff,i);
end

SymbolsFSpreadFilt(S,NrOfSymbols) = single(eps*i);
SymbolsFSpreadOffFilt(S,NrOfSymbols) = single(eps*i);
for i=1:NrOfSymbols
    SymbolsFSpreadFilt(:,i) = ...
        real(SymbolsFSpread(:,i)).*FB_odd +1i*imag(SymbolsFSpread(:,i)).*FB_even;
                                        
    SymbolsFSpreadOffFilt(:,i)= ...
        1i*imag(SymbolsFSpreadOff(:,i)).*FB_odd + real(SymbolsFSpreadOff(:,i)).*FB_even;  
end

%% Time domain
Scale = single(1);
% Oversampling
SymbolsTOv = single(zeros(S * OV,NrOfSymbols));
SymbolsTOv (1:S/2,:)= SymbolsFSpreadFilt(1:S/2,:);
SymbolsTOv (end-S/2+1:end,:)= SymbolsFSpreadFilt(S/2+1:end,:);
SymbolsTOv = ifft(SymbolsTOv) * Scale * OV;

SymbolsTOffOv = single(zeros(S*OV,NrOfSymbols));
SymbolsTOffOv (1:S/2,:)= SymbolsFSpreadOffFilt(1:S/2,:);
SymbolsTOffOv (end-S/2+1:end,:)= SymbolsFSpreadOffFilt(S/2+1:end,:);
SymbolsTOffOv = ifft(SymbolsTOffOv) * Scale* OV;



% signal memory allocation;
signalTx = single(zeros((N/2 + N * NrOfSymbols + (K-1)*N)*OV,1));
for i=1:NrOfSymbols
    
    index_start = 1 + (i-1)*N*OV;
    index_end   =  (i-1) * N*OV + S*OV;
        
    % Normal symbols
    signalTx(index_start : index_end) = signalTx(index_start : index_end)+ SymbolsTOv(:,i);
    % Offseted symbols - complex values
    signalTx(index_start + Offset*OV : index_end + Offset*OV) = ...
        signalTx(index_start + Offset*OV:index_end + Offset*OV) + SymbolsTOffOv(:,i);

    
end

%% Save everything into the return object
Modulated.NrOfSymbols = NrOfSymbols;
Modulated.NrOfExtModSymbs = NrOfExtModSymbs;
Modulated.SymbolsF = SymbolsF;
Modulated.SymbolsFOv = SymbolsFOv;
Modulated.SymbolsFSpread = SymbolsFSpread;
Modulated.SymbolsFOvOff = SymbolsFOvOff;
Modulated.SymbolsFSpreadOff = SymbolsFSpreadOff;
Modulated.SymbolsFSpreadFilt = SymbolsFSpreadFilt;
Modulated.SymbolsFSpreadOffFilt = SymbolsFSpreadOffFilt;
Modulated.SymbolsTOv = SymbolsTOv;
Modulated.SymbolsTOffOv = SymbolsTOffOv;
Modulated.Scale = Scale;
Modulated.signalTx = signalTx;
Modulated.Es = mean(abs(signalTx).^2);
end

