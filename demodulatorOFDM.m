function [ Demodulated ] = demodulatorOFDM(dataSerial,Param)

sizeOfFFT = Param.N;
numOfSym = Param.numOfSym;
overSampling = Param.OV;
cp = Param.CP;
CarrierIndexes = Param.CarrierIndexes;
  
Demodulated.dataSerial = reshape(dataSerial, overSampling*(sizeOfFFT+cp),numOfSym).';
    % 1. Remove cyclic prefix
    if (cp~=0)
     Demodulated.removedCP=Demodulated.dataSerial(:,overSampling*cp+1:end);
    else
     Demodulated.removedCP = Demodulated.dataSerial;
    end
    % 2. Serial to parallel conversion
     Demodulated.dataParallel=Demodulated.removedCP.';
    % 3. Downsampling
     Demodulated.dataParallel = downsample(Demodulated.dataParallel,overSampling);
    % 4. Perform FFT
     Demodulated.downsampledFFTdata = fft(Demodulated.dataParallel);
    % 5. Parallel to serial conversion
     Demodulated.downsampledFFTdataSerial = Demodulated.downsampledFFTdata.';
     
     Demodulated.Symbols = Demodulated.downsampledFFTdataSerial(:,CarrierIndexes);
     
end

