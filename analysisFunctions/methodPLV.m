function [rawPLV, shuffledPLV] = methodPLV (amp, phase)
% Phase locking value, Cohen 2008

% number of Surrogates
nSurrogates = size(amp,2); % 1- filtered raw signal; >1-filtered shuffle signal

if nSurrogates ==1  % PLV for filtered raw signal with one time series
    % Phase of amplitude envelope    
    phaseAmplitudeEnvelope = angle(hilbert(detrend(amp(:)))); 
    rawPLV = abs(mean(exp(1i*(phase(:) - phaseAmplitudeEnvelope(:)))));
    shuffledPLV = 0;
else % PLV for filtered shuffled signal with nSurrogates
    shuffledPLV = zeros(1, nSurrogates);
    parfor nSurr = 1:nSurrogates        
        % Phase of amplitude envelope    
        phaseAmplitudeEnvelope = angle(hilbert(detrend(amp(:, nSurr)))); 
        shuffledPLV(1, nSurr) = abs(mean(exp(1i*(phase(:, nSurr) - phaseAmplitudeEnvelope(:)))));    
    end
    rawPLV = 0;  
end

end
