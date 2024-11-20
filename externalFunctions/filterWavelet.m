function [filteredAmp, filteredPhase] = filterWavelet(amp, ampFreq, phase, phaseFreq, samplingFreq, timeWin)

% Get number of bins as per phase frequencies and amplitude frequencies
phaseBins = ceil((max(phaseFreq) - min(phaseFreq))/(diff(phaseFreq(1:2))));
ampBins = ceil((max(ampFreq) - min(ampFreq))/(diff(ampFreq(1:2))));

% Details of data
numTimes = length(timeWin);
numTrials = size(amp,2);

% Frequency limits 
phaseFreqLimits = [min(phaseFreq) max(phaseFreq)];
ampFreqLimits = [min(ampFreq) max(ampFreq)];
phaseFreqBandwidth = diff(phaseFreq(1:2));
ampFreqBandwidth = diff(ampFreq(1:2));

% Cells to store filtered signals
filteredPhase = cell(1,phaseBins);
filteredAmp = cell(1,ampBins);

% for nPhaseBin=1:phaseBins
%     filteredPhase{1,nPhaseBin} = zeros(numTimes, numTrials);
% end
% 
% for nAmpBin=1:ampBins
%     filteredAmp{1,nAmpBin} = zeros(numTimes, numTrials);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for Filtering
waveletWidth = 7; % width of morlet wavelet

% Filter amplitude signals at different amplitude frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
for nAmpBin = 1:ampBins    
    upperBinFreq = ampFreqLimits(1) + (nAmpBin * ampFreqBandwidth);
    lowerBinFreq = upperBinFreq - ampFreqBandwidth;     
    for nTrials = 1:numTrials
            tmpFilteredSigAll = ampvec((lowerBinFreq + floor((upperBinFreq - lowerBinFreq)/2)), amp(:, nTrials), samplingFreq, waveletWidth);
            tmpFilteredSig(:,nTrials) = tmpFilteredSigAll(:, timeWin);
    end    
    filteredAmp{1,nAmpBin} = tmpFilteredSig(:); % Conatenated trials per each amplitude bin  
end

% Filter phase signals at different phase frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
for nPhaseBin = 1:phaseBins    
    upperBinFreq = phaseFreqLimits(1) + (nPhaseBin * phaseFreqBandwidth);
    lowerBinFreq = upperBinFreq - phaseFreqBandwidth;     
    for nTrials = 1:numTrials
            tmpFilteredSigAll = phasevec((lowerBinFreq + floor((upperBinFreq - lowerBinFreq)/2)), phase(:, nTrials), samplingFreq, waveletWidth);
            tmpFilteredSig(:,nTrials) = tmpFilteredSigAll(:, timeWin);
    end    
    filteredPhase{1,nPhaseBin} = tmpFilteredSig(:); % Concatenated trials per each phase bin 
end

end
