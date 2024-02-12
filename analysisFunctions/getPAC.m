function [rawPac, normalisedPac, shuffledPac, tval, pval, centerAmpFreq, centerPhaseFreq] = getPAC (amp, ampFreq, phase, phaseFreq, samplingFreq, filterName, epochFrames, phFiltOrder, ...
                                    ampFiltOrder, waveletWidth, pacMethod, nSurrogates, alpha)

% Get number of bins as per phase frequencies and amplitude frequencies
phaseBins = ceil((max(phaseFreq) - min(phaseFreq))/(diff(phaseFreq(1:2))));
ampBins = ceil((max(ampFreq) - min(ampFreq))/(diff(ampFreq(1:2))));

% Get the Phase and amplitude time series based on Filtering methods
if (strcmp(filterName, 'wavelet'))
    [filteredAmp, filteredPhase] = filterWavelet(amp, ampFreq, phase, phaseFreq, samplingFreq, waveletWidth);
elseif (strcmp(filterName, 'fir') || strcmp(filterName, 'mp'))
    [filteredAmp, filteredPhase] = filterFir(amp, ampFreq, phase, phaseFreq, samplingFreq, epochFrames, phFiltOrder, ampFiltOrder);
end

% Generate surrogates by shuffling the signals
if nSurrogates ~= 0
    [shuffledAmp, shuffledPhase] = generateSurrogates (filteredAmp, filteredPhase, nSurrogates);

    % PAC on shuffled signals
    [~, shuffledPac, ~, ~] = getPACallMethods (shuffledAmp, ampFreq, shuffledPhase, phaseFreq, pacMethod, nSurrogates);
    disp('End of surrogate analysis');
end

% Different PAC measures on filetered data
[rawPac, ~, centerAmpFreq, centerPhaseFreq] = getPACallMethods (filteredAmp, ampFreq, filteredPhase, phaseFreq, pacMethod, nSurrogates);

%%%%%%%%% Statistical Analysis between Surrogate and Raw PAC values %%%%%%%%%%%%%%

pval = zeros(ampBins, phaseBins); % array for probaility values
tval = zeros(ampBins, phaseBins); % array for t- values
normalisedPac = zeros(ampBins, phaseBins);

% Compute pvalue by estimating a normal distribution from the surrogate data
for nAmpBin = 1:ampBins
    parfor nPhaseBin = 1:phaseBins
        [surrogateMean, surrogateStd] = normfit(squeeze(shuffledPac(nPhaseBin, nAmpBin, :))); % Estimating a normal distribution
        normalisedPac(nAmpBin, nPhaseBin) = (rawPac(nPhaseBin, nAmpBin) - surrogateMean) / surrogateStd; % Normalize the pacval (z-score)
%         pval(nAmpBin, nPhaseBin) = 1-normcdf(abs(normalisedPac(nAmpBin, nPhaseBin)));  % Compute p_value
    end
end

% Compute pvalue using ttest (Surrogate as null distribution)
for nAmpBin = 1:ampBins
    parfor nPhaseBin = 1:phaseBins
        [~, p, ~, stats] = ttest2 (shuffledPac(nPhaseBin, nAmpBin, :), rawPac(nPhaseBin, nAmpBin, :), 'Tail', 'both', 'Alpha', alpha, 'Vartype', 'equal');
        tval(nAmpBin, nPhaseBin) = stats.tstat;
        pval(nAmpBin, nPhaseBin) = p;    
    end
end
disp('End of statistical analysis');
end