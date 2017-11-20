function Modulated = quant_fixpoint_modulatorFBMC(ModulationSymbols, Param)

%% Parameters
S = Param.S;
N = Param.N;
OV = Param.OV;              % Oversampling 
Offset = Param.Offset;      
D = Param.D;
K = Param.K;                % Overlapping factor
CarrierIndexes = Param.CarrierIndexes;

Modulated.NrOfSymbols = ceil(length(ModulationSymbols)/D);
Modulated.NrOfExtModSymbs = mod(D - mod(length(ModulationSymbols),D),D);

%%
SymbolsF = complex(zeros(D, Modulated.NrOfSymbols));
SymbolsF(CarrierIndexes,:) = reshape([ModulationSymbols;zeros(Modulated.NrOfExtModSymbs,1)],D,Modulated.NrOfSymbols);

odd_indexes  = CarrierIndexes(mod(CarrierIndexes,2)==0);
OddIndexesSpread= (odd_indexes-1)*K+1;
even_indexes = CarrierIndexes(mod(CarrierIndexes,2)==1);
EvenIndexesSpread=(even_indexes-1)*K+1;

% Normal symbols
SymbolsFOv = complex(zeros(S,Modulated.NrOfSymbols));
SymbolsFOv(OddIndexesSpread,:)= real(SymbolsF(odd_indexes,:));
SymbolsFOv(EvenIndexesSpread,:)= SymbolsFOv(EvenIndexesSpread,:)+ 1i*real(SymbolsF(even_indexes,:));

SymbolsFSpread = SymbolsFOv*0;
for i=-K:K
    SymbolsFSpread = SymbolsFSpread + circshift(SymbolsFOv,i);
end

% OFFSETSYMBOLS
SymbolsFOvOff=complex(zeros(S,Modulated.NrOfSymbols));
SymbolsFOvOff(OddIndexesSpread,:)= 1i*imag(SymbolsF(odd_indexes,:));
SymbolsFOvOff(EvenIndexesSpread,:)= SymbolsFOvOff(EvenIndexesSpread,:)+ imag(SymbolsF(even_indexes,:));

SymbolsFSpreadOff = SymbolsFOvOff*0;
for i=-K:K
    SymbolsFSpreadOff = SymbolsFSpreadOff + circshift(SymbolsFOvOff,i);
end

SymbolsFSpreadFilt = complex(zeros(S,Modulated.NrOfSymbols));
SymbolsFSpreadOffFilt = complex(zeros(S,Modulated.NrOfSymbols));
for i=1:Modulated.NrOfSymbols
    SymbolsFSpreadFilt(:,i) = real(SymbolsFSpread(:,i)).*Param.FB_odd +...
                                        1i*imag(SymbolsFSpread(:,i)).*Param.FB_even;
    SymbolsFSpreadOffFilt(:,i)= 1i*imag(SymbolsFSpreadOff(:,i)).*Param.FB_odd+...
                                           real(SymbolsFSpreadOff(:,i)).*Param.FB_even;
    
end

Scale = 1/Param.S;
% Oversampling
SymbolsTOv = complex(zeros(S * OV,Modulated.NrOfSymbols));
SymbolsTOv (1:S/2,:)= SymbolsFSpreadFilt(1:S/2,:);
SymbolsTOv (end-S/2+1:end,:)= SymbolsFSpreadFilt(S/2+1:end,:);
for i=1:Modulated.NrOfSymbols
    SymbolsTOv(:,i) = ifft_ct(SymbolsTOv(:,i),S*OV,1) * Scale * Param.OV;
end

SymbolsTOffOv = complex(zeros(S*OV,Modulated.NrOfSymbols));
SymbolsTOffOv (1:S/2,:)= SymbolsFSpreadOffFilt(1:S/2,:);
SymbolsTOffOv (end-S/2+1:end,:)= SymbolsFSpreadOffFilt(S/2+1:end,:);
for i=1:Modulated.NrOfSymbols
    SymbolsTOffOv(:,i) = ifft_ct(SymbolsTOffOv(:,i),S*OV,1) * Scale * Param.OV;
end

Modulated.Scale = Scale;

% signal memory allocation;
signalTx = complex(zeros((N/2 + N * Modulated.NrOfSymbols + (K-1)*N)*OV,1));
for i=1:Modulated.NrOfSymbols
    
    index_start = 1 + (i-1)*N*OV;
    index_end   =  (i-1) * N*OV + S*OV;
        
    % Normal symbols
    signalTx(index_start : index_end) = signalTx(index_start : index_end)+ SymbolsTOv(:,i);
    % Offseted symbols - complex values
    signalTx(index_start + Offset * OV : index_end + Offset * OV) = signalTx(index_start + Offset*OV:index_end + Offset*OV) + SymbolsTOffOv(:,i)  ;

    
end

Modulated.SymbolsF = SymbolsF;
Modulated.SymbolsFOv = SymbolsFOv;
Modulated.SymbolsFSpread = SymbolsFSpread;
Modulated.SymbolsFOvOff = SymbolsFOvOff;
Modulated.SymbolsFSpreadOff = SymbolsFSpreadOff;
Modulated.SymbolsFSpreadFilt = SymbolsFSpreadFilt;
Modulated.SymbolsFSpreadOffFilt = SymbolsFSpreadOffFilt;
Modulated.SymbolsTOv = SymbolsTOv;
Modulated.SymbolsTOffOv = SymbolsTOffOv;
Modulated.signalTx = signalTx;

Modulated.Es = mean(abs(Modulated.signalTx).^2);
end


