function [ Demodulated ] = demodulatorFBMC(RxSignal, Param)

Demodulated.singalRx = RxSignal;

OV = Param.OV;
N = Param.N;
S = Param.S;
D = Param.D;
K = Param.K;
Offset = Param.Offset;
CarrierIndexes = Param.CarrierIndexes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%THETA
k = 1:N;

THETA = exp(-1i*pi/2*k);
THETA2 = exp((-1i*pi/2*(k+1)));
ComplexSymbols(1+(Param.NrOfSymbols-1)*D:Param.NrOfSymbols*D,1) = eps*1i;

for i = 1: Param.NrOfSymbols
        index_start     = 1 + (i-1) * N*OV;
        index_end       = index_start + S*OV -1;
        %%%%%%%%
  
    % Downsampling
    Symbol_t = RxSignal(index_start:OV:index_end);
    Symbol_f =fft(Symbol_t);
    Symbol_f_comp = Symbol_f;
   
    Symbol_offset_t = RxSignal(index_start + Offset*OV : OV : index_end + Offset * OV);
    Symbol_offset_f = fft(Symbol_offset_t);
    Symbol_offset_f_comp = Symbol_offset_f; 
    
   
    % Normal symbol
    data_odd        = real(Symbol_f_comp).*Param.FB_odd;
    data_odd_real   = sum(reshape(circshift(data_odd,0),8,Param.S/8))/4;
    data_odd_real   = data_odd_real.*real(THETA(2:2:N));
    
    data_even       = imag(Symbol_f_comp).*Param.FB_even;
    data_even_real  = sum(reshape(circshift(data_even,4),8,Param.S/8))/4;
    data_even_real  = -1*data_even_real.*imag(THETA(1:2:N));
    
    %Offset symbol
    data_odd        = imag(Symbol_offset_f_comp).*Param.FB_odd;
    data_odd_imag   = sum(reshape(circshift(data_odd,0),8,Param.S/8))/4;
    data_odd_imag   = -1*data_odd_imag.*imag(THETA2(2:2:N));

    data_even       = real(Symbol_offset_f_comp).*Param.FB_even;
    data_even_imag  = sum(reshape(circshift(data_even,4),8,Param.S/8))/4;
    data_even_imag  = data_even_imag.*real(THETA2(1:2:N));
    
    data = zeros(Param.N,1);
    data(2:2:end) = complex(data_odd_real,data_odd_imag);
    data(1:2:end) = complex(data_even_real,data_even_imag);
    
    ComplexSymbols(1+(i-1)*D:i*D,1) = data(CarrierIndexes);
end
    
Demodulated.Symbols = ComplexSymbols;
    
end