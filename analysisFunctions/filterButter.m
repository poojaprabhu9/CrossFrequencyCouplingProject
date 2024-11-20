function [filteredAmp, filteredPhase] = filterButter(amp, ampFreq, phase, phaseFreq, samplingFreq, timeWin)

% Get number of bins as per phase frequencies and amplitude frequencies
phaseBins = ceil((max(phaseFreq) - min(phaseFreq))/(diff(phaseFreq(1:2))));
ampBins = ceil((max(ampFreq) - min(ampFreq))/(diff(ampFreq(1:2))));

% Details of data
% numTimes = length(timeWin);
numTimes = size(amp,1);
numTrials = size(amp,2);

% Frequency bandwidths (or binsize) 
phaseFreqBandwidth = diff(phaseFreq(1:2));
ampFreqBandwidth = diff(ampFreq(1:2));

% Cells to store filtered signals
filteredPhase = cell(1,phaseBins);
filteredAmp = cell(1,ampBins);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter amplitude signals at different amplitude frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
envelopeFilteredSigAll = zeros(numTimes, numTrials);
for nAmpBin = 1:ampBins    
    for nTrials = 1:numTrials
        tmpFilteredSig(:,nTrials) = filterFuncButter(amp(:,nTrials)', samplingFreq, ampFreq(nAmpBin), ampFreq(nAmpBin) + ampFreqBandwidth);
        envelopeFilteredSigAll(:,nTrials) = abs(hilbert(tmpFilteredSig(:,nTrials))); % Magnitude of filtered signal 
    end     
%   to minimise the edge effects, consider only signal corresponding to the analysisPeriod 
    envelopeFilteredSig = envelopeFilteredSigAll(timeWin, :);
    filteredAmp{1,nAmpBin} = envelopeFilteredSig(:); % Concatenate all the trials
end

% Filter phase signals at different phase frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
phaseFilteredSigAll = zeros(numTimes, numTrials);
for nPhaseBin = 1:phaseBins    
    for nTrials = 1:numTrials
        tmpFilteredSig(:,nTrials) = filterFuncButter(phase(:,nTrials)', samplingFreq, phaseFreq(nPhaseBin), phaseFreq(nPhaseBin) + phaseFreqBandwidth);
        phaseFilteredSigAll(:,nTrials) = angle(hilbert(tmpFilteredSig(:,nTrials))); % Phase of filtered signal 
    end    
%      to minimise the edge effects, consider only signal corresponding to the analysisPeriod 
    phaseFilteredSig = phaseFilteredSigAll(timeWin, :);
    filteredPhase{1,nPhaseBin} = phaseFilteredSig(:); % Concatenate all the trials
end

end

function filteredSig = filterFuncButter(sig, samplingFreq, lowCutoff, highCutoff)
    if (lowCutoff >= 1) &&  (lowCutoff <= 100)
        filterOrder = 3;
    else 
        filterOrder = 7;
    end

    normBand = [lowCutoff highCutoff]./(samplingFreq/2);        
%     filterOrder = 3;
    [b,a] = butter(filterOrder,normBand,'bandpass');
    filteredSig = filtfilt(b, a, sig);   

end



