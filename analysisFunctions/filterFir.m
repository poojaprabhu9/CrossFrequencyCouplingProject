function [filteredAmp, filteredPhase, filteredTimeSigAmp, filteredTimeSigPhase] = filterFir(amp, ampFreq, phase, phaseFreq, samplingFreq, timeWin)

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
filteredTimeSigPhase = cell(1,phaseBins);
filteredTimeSigAmp = cell(1,ampBins);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter amplitude signals at different amplitude frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
envelopeFilteredSigAll = zeros(numTimes, numTrials);
for nAmpBin = 1:ampBins    
    for nTrials = 1:numTrials
        tmpFilteredSig(:,nTrials) = filterFun(amp(:,nTrials)', samplingFreq, ampFreq(nAmpBin), ampFreq(nAmpBin) + ampFreqBandwidth);
        envelopeFilteredSigAll(:,nTrials) = abs(hilbert(tmpFilteredSig(:,nTrials))); % Magnitude of filtered signal 
    end 
    filteredTimeSigAmp{1,nAmpBin} = envelopeFilteredSigAll(timeWin, :)';
%   to minimise the edge effects, consider only signal corresponding to the analysisPeriod 
    envelopeFilteredSig = envelopeFilteredSigAll(timeWin, :);
    filteredAmp{1,nAmpBin} = envelopeFilteredSig(:); % Concatenate all the trials
end

% Filter phase signals at different phase frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
phaseFilteredSigAll = zeros(numTimes, numTrials);
for nPhaseBin = 1:phaseBins    
    for nTrials = 1:numTrials
        tmpFilteredSig(:,nTrials) = filterFun(phase(:,nTrials)', samplingFreq, phaseFreq(nPhaseBin), phaseFreq(nPhaseBin) + phaseFreqBandwidth);
        phaseFilteredSigAll(:,nTrials) = angle(hilbert(tmpFilteredSig(:,nTrials))); % Phase of filtered signal 
    end 
    filteredTimeSigPhase{1,nPhaseBin} = phaseFilteredSigAll(timeWin, :)';
    % to minimise the edge effects, consider only signal corresponding to the analysisPeriod 
    phaseFilteredSig = phaseFilteredSigAll(timeWin, :);
    filteredPhase{1,nPhaseBin} = phaseFilteredSig(:); % Concatenate all the trials
end

end

function filteredSig = filterFun(sig, samplingFreq, lowCutoff, highCutoff)
        % Minimum cycles per second
        if (lowCutoff >= 1) &&  (lowCutoff <= 4)
            cyclesPerSec = 1;
        elseif (lowCutoff >= 5) &&  (lowCutoff <= 150) 
            cyclesPerSec = 3; % minimum 3 cycles/sec for lower frequencies
        elseif (lowCutoff >= 150) &&  (lowCutoff <= 500) 
            cyclesPerSec = 10; % minimum 3 cycles/sec for lower frequencies
        end

        % length of the filter in points
        filterOrder = cyclesPerSec *fix(samplingFreq/lowCutoff);

        % 'fir1' is chosen because it gives better attenuation than 'firls'
        filtWeights = fir1(filterOrder, [lowCutoff,  highCutoff]./(samplingFreq/2), 'bandpass');
        filteredSig =  filtfilt(filtWeights, 1, sig);
end



