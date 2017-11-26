function Modulator = quant_single_modulatorFBMC_PPN(ModulationSymbols, Param)
% This polyphase structure does not support oversampling!

%% Parameters
N = single(Param.N);
Offset = Param.Offset;      
D = Param.D;
K = Param.K;
CarrierIndexes = Param.CarrierIndexes;
H_coeffsFB = single(Param.H_coeffsFB);

NrOfSymbols = ceil(length(ModulationSymbols)/D);
NrOfExtModSymbs = mod(D - mod(length(ModulationSymbols),D),D);
SymbolsF(CarrierIndexes,:) = single(reshape([ModulationSymbols;zeros(NrOfExtModSymbs,1)],D,NrOfSymbols));

Theta = single(exp(1i*pi/2*[0:N-1]).');
SymbolsFreal = real(SymbolsF); 
SymbolsFrealtheta = bsxfun(@times,Theta,SymbolsFreal);

Theta = single(exp(1i*pi/2*[1:N]).');
SymbolsFimag = imag(SymbolsF); 
SymbolsFimagtheta = bsxfun(@times,Theta,SymbolsFimag);

%% N - IFFT
SymbolsTreal= single(N*[ifft(SymbolsFrealtheta),zeros(N,K-1)]);
SymbolsTimag= single(N*[ifft(SymbolsFimagtheta),zeros(N,K-1)]);


%% Polyphase structure
pmatrix(N,K) = single(eps*1i);
p0 = ifft(circshift([H_coeffsFB;zeros(K*N-length(H_coeffsFB),1)],-4));
for index=1:N
    pmatrix(index,1:K)=p0([index:N:K*N]);
end

pout_real(N,NrOfSymbols+K-1) = single(eps*1i);
pout_imag(N,NrOfSymbols+K-1) = single(eps*1i);

for index=1:N
pout_real(index,:) = filter(pmatrix(index,:), 1, SymbolsTreal(index,:));
pout_imag(index,:) = filter(pmatrix(index,:), 1, SymbolsTimag(index,:));
end

s_real = pout_real(:).';
s_imag = pout_imag(:).';

signalTx = [s_real(1:end).'; zeros(1,Offset)'];
signalTx(Offset+1:end) = signalTx(Offset+1:end) + s_imag.';

%%Save everything to the return object
Modulator.NrOfSymbols = NrOfSymbols;
Modulator.NrOfExtModSymbs = NrOfExtModSymbs;
Modulator.SymbolsF = SymbolsF;
Modulator.SymbolsFreal = SymbolsFreal;
Modulator.SymbolsFrealtheta = SymbolsFrealtheta;
Modulator.SymbolsFimagtheta = SymbolsFimagtheta;
Modulator.SymbolsTreal = SymbolsTreal;
Modulator.SymbolsTimag = SymbolsTimag;
Modulator.signalTx = signalTx;
Modulator.Es = mean(abs(signalTx).^2);

end

