function [rawPac, shuffledPac, centerAmpFreq, centerPhaseFreq] = getPACallMethods (amp, ampFreq, phase, phaseFreq, pacMethod, nSurrogates)
                                                         
% Number of bins as per given phase and amplitude frequencies
phaseBins = ceil((max(phaseFreq) - min(phaseFreq))/(diff(phaseFreq(1:2))));
ampBins = ceil((max(ampFreq) - min(ampFreq))/(diff(ampFreq(1:2))));

% Center frequencies of each bin
centerPhaseFreq = zeros(phaseBins,1);
centerAmpFreq = zeros(ampBins, 1);

for nPhaseBin = 1:phaseBins
    upperBinFreq = min(phaseFreq) + nPhaseBin*(diff(phaseFreq(1:2)));
    lowerBinFreq = upperBinFreq - (diff(phaseFreq(1:2)));
    centerPhaseFreq(nPhaseBin) =  lowerBinFreq + floor((upperBinFreq - lowerBinFreq)/2);
end

for nAmpBin = 1:ampBins
    upperBinFreq = min(ampFreq) + nAmpBin*(diff(ampFreq(1:2)));
    lowerBinFreq = upperBinFreq - (diff(ampFreq(1:2)));
    centerAmpFreq(nAmpBin) =  lowerBinFreq + floor((upperBinFreq - lowerBinFreq)/2);
end

% initialising the output variables
rawPac = zeros(ampBins, phaseBins);
shuffledPac = zeros(ampBins, phaseBins, nSurrogates);

% Calculate PAC using different methods
for nAmpBin = 1:ampBins
    parfor nPhaseBin = 1:phaseBins
         switch pacMethod
             case 'mvl'
                   [rawPac(nAmpBin, nPhaseBin), shuffledPac(nAmpBin, nPhaseBin, :)] = methodMVL (amp{1, nAmpBin}, phase{1, nPhaseBin});
             case 'plv'
                   [rawPac(nAmpBin, nPhaseBin), shuffledPac(nAmpBin, nPhaseBin, :)] = methodPLV (amp{1, nAmpBin}, phase{1, nPhaseBin});
             case 'glm'
                   [rawPac(nAmpBin, nPhaseBin), shuffledPac(nAmpBin, nPhaseBin, :)] = methodGLM (amp{1, nAmpBin}, phase{1, nPhaseBin});
             case 'klmi'
                   nBins = 18;  
                   [rawPac(nAmpBin, nPhaseBin), shuffledPac(nAmpBin, nPhaseBin, :)] = methodKLMI (amp{1, nAmpBin}, phase{1, nPhaseBin}, nBins);
             case 'nmvl'
                   [rawPac(nAmpBin, nPhaseBin), shuffledPac(nAmpBin, nPhaseBin, :)] = methodnMVL (amp{1, nAmpBin}, phase{1, nPhaseBin});
         end
    end
end

end

