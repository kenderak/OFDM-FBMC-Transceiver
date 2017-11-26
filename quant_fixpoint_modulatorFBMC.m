function Modulated = quant_fixpoint_modulatorFBMC(ModulationSymbols, Param)
%% Fixpoint precision settings
WordLength = 64;
FractionLength = 58;
Fp = fimath('ProductMode','SpecifyPrecision',...
		'ProductWordLength',WordLength,'ProductFractionLength',FractionLength);
Fs = fimath('SumMode','SpecifyPrecision',...
  'SumWordLength',WordLength,'SumFractionLength',FractionLength);

%% Parameters
S = Param.S;
N = Param.N;
OV = Param.OV;           
Offset = Param.Offset;      
D = Param.D;
K = Param.K;                
CarrierIndexes = Param.CarrierIndexes;
FB_odd = fixp(Param.FB_odd);
FB_even = fixp(Param.FB_even);
%% IFFT matrix
W = fixp(complex(zeros(S*OV,S*OV)));
for k = 0:S*OV-1
    for n = 0:S*OV-1
        W(k+1,n+1) = fixp(exp(-1i*2*pi/(S*OV)).^(n*k));
    end
end

Modulated.NrOfSymbols = ceil(length(ModulationSymbols)/D);
Modulated.NrOfExtModSymbs = mod(D - mod(length(ModulationSymbols),D),D);

%%
SymbolsF = fixp(complex(zeros(D, Modulated.NrOfSymbols)));
SymbolsF(CarrierIndexes,:) = reshape([ModulationSymbols;zeros(Modulated.NrOfExtModSymbs,1)],D,Modulated.NrOfSymbols);

odd_indexes  = CarrierIndexes(mod(CarrierIndexes,2)==0);
OddIndexesSpread= (odd_indexes-1)*K+1;
even_indexes = CarrierIndexes(mod(CarrierIndexes,2)==1);
EvenIndexesSpread=(even_indexes-1)*K+1;

% Normal symbols
SymbolsFOv = fixp(complex(zeros(S,Modulated.NrOfSymbols)));
SymbolsFOv(OddIndexesSpread,:)= real(SymbolsF(odd_indexes,:));
SymbolsFOv(EvenIndexesSpread,:)= add(Fs,SymbolsFOv(EvenIndexesSpread,:),mpy(Fp,fixp(1i),real(SymbolsF(even_indexes,:))));

SymbolsFSpread = fixp(complex(zeros(S,Modulated.NrOfSymbols)));
for i=-K:K
    SymbolsFSpread = add(Fs,SymbolsFSpread,circshift(SymbolsFOv,i));
end

% OFFSETSYMBOLS
SymbolsFOvOff=fixp(complex(zeros(S,Modulated.NrOfSymbols)));
SymbolsFOvOff(OddIndexesSpread,:)= mpy(Fp,fixp(1i),imag(SymbolsF(odd_indexes,:)));
SymbolsFOvOff(EvenIndexesSpread,:)= add(Fs,SymbolsFOvOff(EvenIndexesSpread,:),imag(SymbolsF(even_indexes,:)));

SymbolsFSpreadOff = fixp(complex(zeros(S,Modulated.NrOfSymbols)));
for i=-K:K
    SymbolsFSpreadOff = add(Fs,SymbolsFSpreadOff,circshift(SymbolsFOvOff,i));
end

SymbolsFSpreadFilt = fixp(complex(zeros(S,Modulated.NrOfSymbols)));
SymbolsFSpreadOffFilt = fixp(complex(zeros(S,Modulated.NrOfSymbols)));
for i=1:Modulated.NrOfSymbols
    SymbolsFSpreadFilt(:,i) = add(Fs,mpy(Fp,real(SymbolsFSpread(:,i)),FB_odd),...
                                        mpy(Fp,mpy(Fp,fixp(1i),imag(SymbolsFSpread(:,i))),FB_even));
    SymbolsFSpreadOffFilt(:,i)= add(Fs,mpy(Fp,mpy(Fp,fixp(1i),imag(SymbolsFSpreadOff(:,i))),FB_odd),...
                                           mpy(Fp,real(SymbolsFSpreadOff(:,i)),FB_even));
    
end

Scale = fixp(OV/S);
% Oversampling
SymbolsTOv = fixp(complex(zeros(S * OV,Modulated.NrOfSymbols)));
SymbolsTOv (1:S/2,:)= SymbolsFSpreadFilt(1:S/2,:);
SymbolsTOv (end-S/2+1:end,:)= SymbolsFSpreadFilt(S/2+1:end,:);
for i=1:Modulated.NrOfSymbols
    SymbolsTOv(:,i) = fixp(W*SymbolsTOv(:,i)); %ifft_matrix(SymbolsTOv(:,i),W);
end
SymbolsTOv = mpy(Fp,SymbolsTOv, Scale);

SymbolsTOffOv = fixp(complex(zeros(S * OV,Modulated.NrOfSymbols)));
SymbolsTOffOv (1:S/2,:)= SymbolsFSpreadOffFilt(1:S/2,:);
SymbolsTOffOv (end-S/2+1:end,:)= SymbolsFSpreadOffFilt(S/2+1:end,:);
for i=1:Modulated.NrOfSymbols
    SymbolsTOffOv(:,i) = fixp(W*SymbolsTOffOv(:,i)); %ifft_matrix(SymbolsTOffOv(:,i),W);
end
SymbolsTOffOv = mpy(Fp,SymbolsTOffOv, Scale);

% signal memory allocation;
signalTx = fixp(complex(zeros((N/2 + N * Modulated.NrOfSymbols + (K-1)*N)*OV,1)));
for i=1:Modulated.NrOfSymbols
    
    index_start = 1 + (i-1)*N*OV;
    index_end   =  (i-1) * N*OV + S*OV;
        
    % Normal symbols
    signalTx(index_start : index_end) = add(Fs,signalTx(index_start : index_end),SymbolsTOv(:,i));
    % Offseted symbols - complex values
    signalTx(index_start + Offset * OV : index_end + Offset * OV) =...
        add(Fs,signalTx(index_start + Offset*OV:index_end + Offset*OV),SymbolsTOffOv(:,i));

end

%% Save everything into the return object
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
Modulated.Scale = Scale;

Modulated.Es = mean(abs(Modulated.signalTx).^2);
end


