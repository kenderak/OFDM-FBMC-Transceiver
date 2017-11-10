function Modulator = modulatorFBMC(ModulationSymbols, Param)

%% Parameters
S = Param.S;
N = Param.N;
OV = Param.OV; 
Offset = Param.Offset;      
D = Param.D;
K = Param.K;                
CarrierIndexes = Param.CarrierIndexes;

Modulator.NrOfSymbols = ceil(length(ModulationSymbols)/D);
Modulator.NrOfExtModSymbs = mod(D - mod(length(ModulationSymbols),D),D);

%%
Modulator.SymbolsF(CarrierIndexes,:) = reshape([ModulationSymbols;zeros(Modulator.NrOfExtModSymbs,1)],D,Modulator.NrOfSymbols);

odd_indexes  = CarrierIndexes(mod(CarrierIndexes,2)==0);
OddIndexesSpread= (odd_indexes-1)*K+1;
even_indexes = CarrierIndexes(mod(CarrierIndexes,2)==1);
EvenIndexesSpread=(even_indexes-1)*K+1;

%%THETA
k = 1:N;
THETA = exp(1i*pi/2*k).';
THETA2 = exp(1i*pi/2*(k+1)).';

% Normal symbols
Modulator.SymbolsFOv(S,Modulator.NrOfSymbols) = eps*1i;

Modulator.SymbolsFOv(OddIndexesSpread,:)= bsxfun(@times,THETA(odd_indexes),real(Modulator.SymbolsF(odd_indexes,:)));
Modulator.SymbolsFOv(EvenIndexesSpread,:)= Modulator.SymbolsFOv(EvenIndexesSpread,:)+...
    bsxfun(@times,THETA(even_indexes),real(Modulator.SymbolsF(even_indexes,:)));

% Modulator.SymbolsFOv(OddIndexesSpread,:)= real(Modulator.SymbolsF(odd_indexes,:));
% Modulator.SymbolsFOv(EvenIndexesSpread,:)= Modulator.SymbolsFOv(EvenIndexesSpread,:)+...
%     1i*real(Modulator.SymbolsF(even_indexes,:));


Modulator.SymbolsFSpread = Modulator.SymbolsFOv*0;
for i=-K:K
    Modulator.SymbolsFSpread = Modulator.SymbolsFSpread + circshift(Modulator.SymbolsFOv,i);
end

% OFFSETSYMBOLS
Modulator.SymbolsFOvOff(S,Modulator.NrOfSymbols)=eps*1i;

Modulator.SymbolsFOvOff(OddIndexesSpread,:)= bsxfun(@times,THETA2(odd_indexes),imag(Modulator.SymbolsF(odd_indexes,:)));
Modulator.SymbolsFOvOff(EvenIndexesSpread,:)= Modulator.SymbolsFOvOff(EvenIndexesSpread,:)+...
   bsxfun(@times,THETA2(even_indexes),imag(Modulator.SymbolsF(even_indexes,:)));

% Modulator.SymbolsFOvOff(OddIndexesSpread,:)= 1i*imag(Modulator.SymbolsF(odd_indexes,:));
% Modulator.SymbolsFOvOff(EvenIndexesSpread,:)= Modulator.SymbolsFOvOff(EvenIndexesSpread,:)+...
%    imag(Modulator.SymbolsF(even_indexes,:));


Modulator.SymbolsFSpreadOff = Modulator.SymbolsFOvOff*0;
for i=-K:K
    Modulator.SymbolsFSpreadOff = Modulator.SymbolsFSpreadOff + circshift(Modulator.SymbolsFOvOff,i);
end

for i=1:Modulator.NrOfSymbols
    Modulator.SymbolsFSpreadFilt(:,i) = real(Modulator.SymbolsFSpread(:,i)).*Param.FB_odd +...
                                        1i*imag(Modulator.SymbolsFSpread(:,i)).*Param.FB_even;
    Modulator.SymbolsFSpreadOffFilt(:,i)= 1i*imag(Modulator.SymbolsFSpreadOff(:,i)).*Param.FB_odd+...
                                           real(Modulator.SymbolsFSpreadOff(:,i)).*Param.FB_even;
    
end

% Time domain
Scale = 1;

% Oversampling
Modulator.SymbolsTOv = zeros(S * OV,Modulator.NrOfSymbols);
Modulator.SymbolsTOv (1:S/2,:)= Modulator.SymbolsFSpreadFilt(1:S/2,:);
Modulator.SymbolsTOv (end-S/2+1:end,:)= Modulator.SymbolsFSpreadFilt(S/2+1:end,:);
Modulator.SymbolsTOv = ifft(Modulator.SymbolsTOv) * Scale * Param.OV;

Modulator.SymbolsTOffOv = zeros(S*OV,Modulator.NrOfSymbols);
Modulator.SymbolsTOffOv (1:S/2,:)= Modulator.SymbolsFSpreadOffFilt(1:S/2,:);
Modulator.SymbolsTOffOv (end-S/2+1:end,:)= Modulator.SymbolsFSpreadOffFilt(S/2+1:end,:);
Modulator.SymbolsTOffOv = ifft(Modulator.SymbolsTOffOv) * Scale*Param.OV;

Modulator.Scale = Scale;

% signal memory allocation;
Modulator.signalTx = zeros((N/2 + N * Modulator.NrOfSymbols + (3)*N)*OV,1);
for i=1:Modulator.NrOfSymbols
    
    index_start = 1 + (i-1)*N*OV;
    index_end   =  (i-1) * N*OV + S*OV;
        
    % Normal symbols
    Modulator.signalTx(index_start : index_end) = Modulator.signalTx(index_start : index_end)+ Modulator.SymbolsTOv(:,i);
    % Offseted symbols - complex values
    Modulator.signalTx(index_start + Offset * OV : index_end + Offset * OV) = Modulator.signalTx(index_start + Offset*OV:index_end + Offset*OV) + Modulator.SymbolsTOffOv(:,i)  ;

    
end
Modulator.Es = mean(abs(Modulator.signalTx).^2);

