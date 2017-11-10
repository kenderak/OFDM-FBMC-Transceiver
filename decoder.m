function [ decodedData ] = decoder( demappedData, codingTechnique )
switch codingTechnique
    case 'None'
        decodedData = demappedData;

end

