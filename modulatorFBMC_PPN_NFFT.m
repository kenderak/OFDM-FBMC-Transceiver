function Modulator = modulatorFBMC_PPN_NFFT(ModulationSymbols, Param)

%% Parameters
S = Param.S;
N = Param.N;
OV = Param.OV;              % Oversampling 
Offset = Param.Offset;      
D = Param.D;
K = Param.K;                % Overlapping factor
CarrierIndexes = Param.CarrierIndexes;

Modulator.NrOfSymbols = ceil(length(ModulationSymbols)/D);
Modulator.NrOfExtModSymbs = mod(D - mod(length(ModulationSymbols),D),D);

%%
Modulator.SymbolsF(CarrierIndexes,:) = reshape([ModulationSymbols;zeros(Modulator.NrOfExtModSymbs,1)],D,Modulator.NrOfSymbols);

%% N - IFFT
Modulator.SymbolsT = N*[ifft(Modulator.SymbolsF), zeros(N,3)];
%% Signal separation
[H, G] = separate(Modulator.SymbolsT);
%% Circular timeshift by N/4
Modulator.H = circshift(H,[-N/4 0]);
Modulator.G = circshift(1i*G,[-N/4 0]);

%% Filter
pmatrix(N,K) = eps*1i;
p0 = ifft(circshift([Param.H_coeffsFB;zeros(K*N-length(Param.H_coeffsFB),1)],-4));
for index=1:N
pmatrix(index,1:K)=p0([index:N:K*N]);
end

%% Polyphase structure
pout_real(N,Modulator.NrOfSymbols+K-1) = eps*1i;
pout_imag(N,Modulator.NrOfSymbols+K-1) = eps*1i;

for index=1:N
pout_real(index,:) = filter(pmatrix(index,:), 1, Modulator.H(index,:));
pout_imag(index,:) = filter(pmatrix(index,:), 1, Modulator.G(index,:));
end

s_real=[pout_real(:,1).'];
s_imag=[pout_imag(:,1).'];
for index=2:Modulator.NrOfSymbols+3
    s_real  = [s_real, pout_real(:,index).']; 
    s_imag  = [s_imag, pout_imag(:,index).']; 
end

Modulator.signalTx = [s_real(1:end).'; zeros(1,N/2)'];
Modulator.signalTx(N/2+1:end) = Modulator.signalTx(N/2+1:end) + s_imag.';

Modulator.Es = mean(abs(Modulator.signalTx).^2);
end

