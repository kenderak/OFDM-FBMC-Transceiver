function [ Demodulated ] = demodulatorFBMC(RxSignal, Param)

Demodulated.Rx = RxSignal;

OV = Param.OV;
N = Param.N;
S = Param.S;
D = Param.D;
Offset = Param.Offset;
CarrierIndexes = Param.CarrierIndexes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
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

    data_even       = imag(Symbol_f_comp).*Param.FB_even;
    data_even_real  = sum(reshape(circshift(data_even,4),8,Param.S/8))/4;
    
    %Offset symbol
    data_odd       = imag(Symbol_offset_f_comp).*Param.FB_odd;
    data_odd_imag   = sum(reshape(circshift(data_odd,0),8,Param.S/8))/4;

    data_even       = real(Symbol_offset_f_comp).*Param.FB_even;
    data_even_imag  = sum(reshape(circshift(data_even,4),8,Param.S/8))/4;
    
    data = zeros(Param.N,1);
    data(2:2:end) = complex(data_odd_real,data_odd_imag);
    data(1:2:end) = complex(data_even_real,data_even_imag);
    
    ComplexSymbols(1+(i-1)*D:i*D,1) = data(CarrierIndexes);
end
    
Demodulated.Symbols = ComplexSymbols;
    
end