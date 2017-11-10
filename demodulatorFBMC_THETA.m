function [ Demodulated ] = demodulatorFBMC_THETA(RxSignal, Param)

Demodulated.Rx = RxSignal;

OV = Param.OV;
N = Param.N;
S = Param.S;
D = Param.D;
Offset = Param.Offset;
CarrierIndexes = Param.CarrierIndexes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

ComplexSymbolsF(D,Param.NrOfSymbols)=1i*eps; 

for i = 1: Param.NrOfSymbols
    
    index_start     = 1 + (i-1) * N*OV;
    index_end       = index_start + S*OV -1;
 
  
    % Downsampling
    Symbol_t = RxSignal(index_start:OV:index_end);
    Symbol_f =fft(Symbol_t);
    Symbol_f_comp = Symbol_f;
   
    Symbol_offset_t = RxSignal(index_start + Offset*OV : OV : index_end + Offset * OV);
    Symbol_offset_f = fft(Symbol_offset_t);
    Symbol_offset_f_comp = Symbol_offset_f; 
    
   
    % Normal symbol
    data_odd        = imag(Symbol_f_comp).*Param.FB_odd;
    data_odd_real   = sum(reshape(circshift(data_odd,0),8,Param.S/8))/4;

    data_even       = real(Symbol_f_comp).*Param.FB_even;
    data_even_real  = sum(reshape(circshift(data_even,4),8,Param.S/8))/4;
    
    %Offset symbol
    data_odd        = real(Symbol_offset_f_comp).*Param.FB_odd;
    data_odd_imag   = sum(reshape(circshift(data_odd,0),8,Param.S/8))/4;

    data_even       = imag(Symbol_offset_f_comp).*Param.FB_even;
    data_even_imag  = sum(reshape(circshift(data_even,4),8,Param.S/8))/4;
    
    %% THETA PART
    data = zeros(Param.N,1);
    Theta = exp(1i*pi/2*[0:N-1]);
    Thetakp1 = exp(1i*pi/2*[1:N]);
    data_odd_realtheta = real(data_odd_real.*Theta(1:2:end));
    data_odd_imagtheta = real(data_odd_imag.*Thetakp1(2:2:end));
    
    data_even_realtheta = real(data_even_real.*imag(Theta(2:2:end)));
    data_even_imagtheta = real(data_even_imag.*imag(Thetakp1(1:2:end)));
    
    %%
    data(2:2:end) = complex(data_odd_realtheta,data_odd_imagtheta);
    data(1:2:end) = complex(data_even_realtheta,data_even_imagtheta);
    
    ComplexSymbolsF(:,i) = data(CarrierIndexes);
    
end

Demodulated.ComplexSymbolsF = ComplexSymbolsF;    
Demodulated.Symbols = ComplexSymbolsF(:);
    
end