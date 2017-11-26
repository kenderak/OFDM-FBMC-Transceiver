function [converted] = fixp(ValueToConvert)
%% Fixpoint precision
WordLength = 64;
FractionLength = 58;
Sign = 1;
    
converted = fi(ValueToConvert, Sign, WordLength, FractionLength);

end