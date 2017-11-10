function Modulated = modulatorFBMC(ModulationSymbols, Param)

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
Modulated.SymbolsF(CarrierIndexes,:) = reshape([ModulationSymbols;zeros(Modulated.NrOfExtModSymbs,1)],D,Modulated.NrOfSymbols);

odd_indexes  = CarrierIndexes(mod(CarrierIndexes,2)==0);
OddIndexesSpread= (odd_indexes-1)*K+1;
even_indexes = CarrierIndexes(mod(CarrierIndexes,2)==1);
EvenIndexesSpread=(even_indexes-1)*K+1;

% Normal symbols
Modulated.SymbolsFOv(S,Modulated.NrOfSymbols) = eps*1i;
Modulated.SymbolsFOv(OddIndexesSpread,:)= real(Modulated.SymbolsF(odd_indexes,:));
Modulated.SymbolsFOv(EvenIndexesSpread,:)= Modulated.SymbolsFOv(EvenIndexesSpread,:)+ 1i*real(Modulated.SymbolsF(even_indexes,:));

Modulated.SymbolsFSpread = Modulated.SymbolsFOv*0;
for i=-K:K
    Modulated.SymbolsFSpread = Modulated.SymbolsFSpread + circshift(Modulated.SymbolsFOv,i);
end

% OFFSETSYMBOLS
Modulated.SymbolsFOvOff(S,Modulated.NrOfSymbols)=eps*1i;
Modulated.SymbolsFOvOff(OddIndexesSpread,:)= 1i*imag(Modulated.SymbolsF(odd_indexes,:));
Modulated.SymbolsFOvOff(EvenIndexesSpread,:)= Modulated.SymbolsFOvOff(EvenIndexesSpread,:)+ imag(Modulated.SymbolsF(even_indexes,:));

Modulated.SymbolsFSpreadOff = Modulated.SymbolsFOvOff*0;
for i=-K:K
    Modulated.SymbolsFSpreadOff = Modulated.SymbolsFSpreadOff + circshift(Modulated.SymbolsFOvOff,i);
end

for i=1:Modulated.NrOfSymbols
    Modulated.SymbolsFSpreadFilt(:,i) = real(Modulated.SymbolsFSpread(:,i)).*Param.FB_odd +...
                                        1i*imag(Modulated.SymbolsFSpread(:,i)).*Param.FB_even;
    Modulated.SymbolsFSpreadOffFilt(:,i)= 1i*imag(Modulated.SymbolsFSpreadOff(:,i)).*Param.FB_odd+...
                                           real(Modulated.SymbolsFSpreadOff(:,i)).*Param.FB_even;
    
end

% Time domain
Scale = 1;
Modulated.SymbolsT      = ifft(Modulated.SymbolsFSpreadFilt)*Scale;
Modulated.SymbolsTOff   = ifft(Modulated.SymbolsFSpreadOffFilt)*Scale;

% Oversampling
Modulated.SymbolsTOv = zeros(S * OV,Modulated.NrOfSymbols);
Modulated.SymbolsTOv (1:S/2,:)= Modulated.SymbolsFSpreadFilt(1:S/2,:);
Modulated.SymbolsTOv (end-S/2+1:end,:)= Modulated.SymbolsFSpreadFilt(S/2+1:end,:);
Modulated.SymbolsTOv = ifft(Modulated.SymbolsTOv) * Scale * Param.OV;

Modulated.SymbolsTOffOv = zeros(S*OV,Modulated.NrOfSymbols);
Modulated.SymbolsTOffOv (1:S/2,:)= Modulated.SymbolsFSpreadOffFilt(1:S/2,:);
Modulated.SymbolsTOffOv (end-S/2+1:end,:)= Modulated.SymbolsFSpreadOffFilt(S/2+1:end,:);
Modulated.SymbolsTOffOv = ifft(Modulated.SymbolsTOffOv) * Scale*Param.OV;

Modulated.Scale = Scale;

% signal memory allocation;
Modulated.signalTx = zeros((N/2 + N * Modulated.NrOfSymbols + (K-1)*N)*OV,1);
for i=1:Modulated.NrOfSymbols
    
    index_start = 1 + (i-1)*N*OV;
    index_end   =  (i-1) * N*OV + S*OV;
        
    % Normal symbols
    Modulated.signalTx(index_start : index_end) = Modulated.signalTx(index_start : index_end)+ Modulated.SymbolsTOv(:,i);
    % Offseted symbols - complex values
    Modulated.signalTx(index_start + Offset * OV : index_end + Offset * OV) = Modulated.signalTx(index_start + Offset*OV:index_end + Offset*OV) + Modulated.SymbolsTOffOv(:,i)  ;

    
end
Modulated.Es = mean(abs(Modulated.signalTx).^2);
end

