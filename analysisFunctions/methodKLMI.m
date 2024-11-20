
function [rawKLMI,angleKLMI,shuffledKLMI, meanAmp] = methodKLMI (amp, phase, nBins)
% Kullback–Leibler Modulation Index, Tort et PLVal., 2008 and 2010

% this variable will get the beginning (not the center) of each phase bin (in rads)
position=zeros(1,nBins); 
winsize = 2*pi/nBins;
for nBin=1:nBins 
    position(nBin) = -pi+(nBin-1) * winsize; 
end

% number of Surrogates
nSurrogates = size(amp,2); % 1- filtered raw signal; >1-filtered shuffle signal

if nSurrogates ==1  % KLMI for filtered raw signal with one time series
    meanAmp=zeros(1,nBins); 
    for nBin=1:nBins   
        I = find(phase(:) <  position(nBin) + winsize & phase(:) >=  position(nBin));
        Amp = amp(:);
        meanAmp(nBin) = mean(Amp(I)); 
    end      
    % Quantifying the amount of amp modulation using normalized entropy index (Tort et al PNAS 2008)    
    rawKLMI = (log(nBins) - (-sum((meanAmp/sum(meanAmp)) .* log((meanAmp/sum(meanAmp)))))) / log(nBins);

    % peakangle
    normMeanAmp = meanAmp/(sum(meanAmp)); % normalise
    [~, idx] = max(normMeanAmp);
    angleKLMI = idx*winsize - (winsize/2);
    shuffledKLMI = 0;
else % KLMI for filtered shuffled signal with nSurrogates
    shuffledKLMI = zeros(1, nSurrogates);
    parfor nSurr = 1:nSurrogates        
        % Compute the mean amplitude in each phase
        meanAmp=zeros(1,nBins); 
        for nBin=1:nBins 
            I = find(phase(:, nSurr) <  position(nBin) + winsize & phase(:, nSurr) >=  position(nBin));
            Amp = amp(:, nSurr);
            meanAmp(nBin) = mean(Amp(I)); 
        end 
        % Quantifying the amount of amp modulation using normalized entropy index (Tort et al PNAS 2008)    
        shuffledKLMI(1, nSurr) = (log(nBins) - (-sum((meanAmp/sum(meanAmp)) .* log((meanAmp/sum(meanAmp))))))/log(nBins);

         % peakangle
        normMeanAmp = meanAmp/(sum(meanAmp)); % normalise
        [~, idx] = max(normMeanAmp);
        shuffledAngleKLMI(1, nSurr) = idx*winsize - (winsize/2);
    end
    rawKLMI = 0;  
end

end

