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

Scale = 1/Param.S;
% Oversampling
Modulated.SymbolsTOv = zeros(S * OV,Modulated.NrOfSymbols);
Modulated.SymbolsTOv (1:S/2,:)= Modulated.SymbolsFSpreadFilt(1:S/2,:);
Modulated.SymbolsTOv (end-S/2+1:end,:)= Modulated.SymbolsFSpreadFilt(S/2+1:end,:);
for i=1:Modulated.NrOfSymbols
    Modulated.SymbolsTOv(:,i) = ifft_ct(Modulated.SymbolsTOv(:,i),S*OV,1) * Scale * Param.OV;
end

Modulated.SymbolsTOffOv = zeros(S*OV,Modulated.NrOfSymbols);
Modulated.SymbolsTOffOv (1:S/2,:)= Modulated.SymbolsFSpreadOffFilt(1:S/2,:);
Modulated.SymbolsTOffOv (end-S/2+1:end,:)= Modulated.SymbolsFSpreadOffFilt(S/2+1:end,:);
for i=1:Modulated.NrOfSymbols
    Modulated.SymbolsTOffOv(:,i) = ifft_ct(Modulated.SymbolsTOffOv(:,i),S*OV,1) * Scale * Param.OV;
end

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

%% Cooley-Tukey IFFT algorithm
function Xk= ifft_ct(xn,N,s)
    if (N <= 1)
        Xk(1) = xn(1);
        return
    else
        % Evens!
        xn_even = xn(1:2:N);
        Xk(1:N/2)   = ifft_ct(xn_even,N/2,2*s);
        
        % Odds!
        xn_odd = xn(2:2:N);
        Xk(N/2+1:N) = ifft_ct(xn_odd,N/2,2*s);
        
        for k = 0:(N/2)-1
            t = Xk(k+1);
            Xk(k+1) = t + exp(2*pi*1i*k/N) * Xk(k+1+N/2);
            Xk(k+1+N/2) = t - exp(2*pi*1i*k/N) * Xk(k+1+N/2);
        end
    end
end

end


