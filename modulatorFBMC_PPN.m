function Modulator = modulatorFBMC_PPN(ModulationSymbols, Param)

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

%%THETA
k = 0:N-1;
THETA = exp(1i*pi/2*k).';
THETA2 = exp(1i*pi/2*(k+1)).';

%%
Modulator.SymbolsF(CarrierIndexes,:) = reshape([ModulationSymbols;zeros(Modulator.NrOfExtModSymbs,1)],D,Modulator.NrOfSymbols);

odd_indexes  = CarrierIndexes(mod(CarrierIndexes,2)==0);
OddIndexesSpread= (odd_indexes-1)+1;
even_indexes = CarrierIndexes(mod(CarrierIndexes,2)==1);
EvenIndexesSpread=(even_indexes-1)+1;


Modulator.SymbolsFreal(N,Modulator.NrOfSymbols) = eps*1i;
Modulator.SymbolsFimag(N,Modulator.NrOfSymbols) = eps*1i;

Modulator.SymbolsFreal = zeros(N,Modulator.NrOfSymbols);
Modulator.SymbolsFreal = real(Modulator.SymbolsF);
Modulator.SymbolsFreal(EvenIndexesSpread,:) = Modulator.SymbolsFreal(EvenIndexesSpread,:); % *1i if theta not exist

Modulator.SymbolsFimag = zeros(N,Modulator.NrOfSymbols);
Modulator.SymbolsFimag = imag(Modulator.SymbolsF);
Modulator.SymbolsFimag(OddIndexesSpread,:) = Modulator.SymbolsFimag(OddIndexesSpread,:); % *1i if theta2 not exist

% Multiply by Theta k and k+1
Modulator.SymbolsFreal = bsxfun(@times,THETA,Modulator.SymbolsFreal);
Modulator.SymbolsFimag = bsxfun(@times,THETA2,Modulator.SymbolsFimag);

%% N - IFFT
Modulator.SymbolsTreal= N*[ifft(Modulator.SymbolsFreal),zeros(N,3)];
Modulator.SymbolsTimag= N*[ifft(Modulator.SymbolsFimag),zeros(N,3)];

pmatrix(N,K) = eps*1i;
p0 = ifft(circshift([Param.H_coeffsFB;zeros(K*N-length(Param.H_coeffsFB),1)],-4));
for index=1:N
pmatrix(index,1:K)=p0([index:N:K*N]);
end

%% Polyphase structure
pout_real(N,Modulator.NrOfSymbols+K-1) = eps*1i;
pout_imag(N,Modulator.NrOfSymbols+K-1) = eps*1i;

for index=1:N
pout_real(index,:) = filter(pmatrix(index,:), 1, Modulator.SymbolsTreal(index,:));
pout_imag(index,:) = filter(pmatrix(index,:), 1, Modulator.SymbolsTimag(index,:));
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
