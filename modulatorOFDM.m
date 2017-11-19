function [ Modulated ] = modulatorOFDM(ModulationSymbols,Param)

sizeOfFFT = Param.N;
numOfCarrier = Param.D;
numOfSym = Param.numOfSym;
overSampling = Param.OV;
cp = Param.CP;
CarrierIndexes = Param.CarrierIndexes;

% Create the mapped symbols into a matrix nad insert pilots
Modulated.dataSerial(:,CarrierIndexes) = reshape(ModulationSymbols, numOfSym, numOfCarrier );

 % 1. Serial to parallel conversion
 Modulated.dataParallel=Modulated.dataSerial.';
 % 2. Oversampling with interpolation
 if (overSampling > 1)
    Modulated.OVdataParallel = overSampling*[Modulated.dataParallel(1,:); ...
        Modulated.dataParallel(2:sizeOfFFT/2+1,:); zeros(sizeOfFFT*overSampling-sizeOfFFT,numOfSym);...
        Modulated.dataParallel(sizeOfFFT/2+2:end,:)];
 else
    Modulated.OVdataParallel = overSampling*Modulated.dataParallel;
 end
 % 3. Perform IFFT
 Modulated.IFFT_data=ifft(Modulated.OVdataParallel);
 % 4. Parallel to serial conversion
 Modulated.IFFT_data_serial=Modulated.IFFT_data.';
 % 5. Add cyclic prefix
 if (cp~=0)
    cyclicPrefix=Modulated.IFFT_data_serial(:,end-overSampling*cp+1:end);
    OFDMmodulated = [cyclicPrefix Modulated.IFFT_data_serial].';
 else
    OFDMmodulated = Modulated.IFFT_data_serial.';
 end
 Modulated.signalTx = OFDMmodulated(:);
 Modulated.Es = mean(abs(Modulated.signalTx).^2);
end

