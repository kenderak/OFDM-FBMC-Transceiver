function Modulator = modulatorFBMC_PPN_NFFT(ModulationSymbols, Param)

%% Parameters
N               = Param.N;
Offset          = Param.Offset;      
D               = Param.D;
K               = Param.K;
CarrierIndexes  = Param.CarrierIndexes;

NrOfSymbols     = ceil(length(ModulationSymbols)/D);
NrOfExtModSymbs = mod(D - mod(length(ModulationSymbols),D),D);

%%
SymbolsF(CarrierIndexes,:) = reshape([ModulationSymbols;zeros(NrOfExtModSymbs,1)],D,NrOfSymbols);

%% N - IFFT
SymbolsT = N*[ifft(SymbolsF), zeros(N,K-1)];

%% Signal separation
[H, G] = separate(SymbolsT);

%% Circular timeshift by N/4
H = circshift(H,[-N/4 0]);
G = circshift(1i*G,[-N/4 0]);

%% Filter for polyphase network
pmatrix(N,K) = eps*1i;
p0 = ifft(circshift([Param.H_coeffsFB;zeros(K*N-length(Param.H_coeffsFB),1)],-4));
for index=1:N
pmatrix(index,1:K)=p0([index:N:K*N]);
end

%% Polyphase structure
pout_real(N,NrOfSymbols+K-1) = eps*1i;
pout_imag(N,NrOfSymbols+K-1) = eps*1i;

for index=1:N
pout_real(index,:) = filter(pmatrix(index,:), 1, H(index,:));
pout_imag(index,:) = filter(pmatrix(index,:), 1, G(index,:));
end

s_real = pout_real(:).';
s_imag = pout_imag(:).';

%% FBMC signal generating
signalTx = [s_real(1:end).'; zeros(1,Offset)'];
signalTx(Offset+1:end) = signalTx(Offset+1:end) + s_imag.';

%% Save everything for the return object
Modulator.NrOfSymbols = NrOfSymbols;
Modulator.NrOfExtModSymbs = NrOfExtModSymbs;
Modulator.SymbolsF = SymbolsF;
Modulator.SymbolsT = SymbolsT;
Modulator.H = H;
Modulator.G = G;
Modulator.signalTx = signalTx;
Modulator.Es = mean(abs(Modulator.signalTx).^2);
end

