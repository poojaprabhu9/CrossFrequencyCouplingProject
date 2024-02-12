function [filteredAmp, filteredPhase] = filterFir(amp, ampFreq, phase, phaseFreq, samplingFreq, epochFrames, phFiltOrder, ampFiltOrder)

% Get number of bins as per phase frequencies and amplitude frequencies
phaseBins = ceil((max(phaseFreq) - min(phaseFreq))/(diff(phaseFreq(1:2))));
ampBins = ceil((max(ampFreq) - min(ampFreq))/(diff(ampFreq(1:2))));

% Details of data
numTimes = size(amp,1);
numTrials = size(amp,2);

% Frequency bandwidths (or binsize) 
phaseFreqBandwidth = diff(phaseFreq(1:2));
ampFreqBandwidth = diff(ampFreq(1:2));

% Cells to store filtered signals
filteredPhase = cell(1,phaseBins);
filteredAmp = cell(1,ampBins);

% Additional filter parameters
revFilt = 0; % bandpass to noth filtering
firType = 'fir1'; % for larger attenuation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter amplitude signals at different amplitude frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
for nAmpBin = 1:ampBins    
    for nTrials = 1:numTrials
        tmpFilteredSig(:,nTrials) = eegfilt(amp(:,nTrials)', samplingFreq, ampFreq(nAmpBin), ampFreq(nAmpBin) + ampFreqBandwidth, epochFrames, ampFiltOrder, revFilt, firType);
    end    
    filteredAmp{1,nAmpBin} = abs(hilbert(tmpFilteredSig(:))); % Magnitude of filtered signal having concatenated trials in each amplitude bin  
end

% Filter phase signals at different phase frequencies
tmpFilteredSig = zeros(numTimes, numTrials);
for nPhaseBin = 1:phaseBins    
    for nTrials = 1:numTrials
        tmpFilteredSig(:,nTrials) = eegfilt(phase(:,nTrials)', samplingFreq, phaseFreq(nPhaseBin), phaseFreq(nPhaseBin) + phaseFreqBandwidth, epochFrames, phFiltOrder, revFilt, firType);
    end    
    filteredPhase{1,nPhaseBin} = angle(hilbert(tmpFilteredSig(:))); % Phase of filtered signal having concatenated trials in each amplitude bin   
end

end


