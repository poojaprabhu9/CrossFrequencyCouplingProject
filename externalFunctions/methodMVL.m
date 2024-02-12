function [rawMVL, shuffledMVL] = methodMVL (amp, phase)
% Mean-Vector length modulation index, Canolty et al., 2006

% number of Surrogates
nSurrogates = size(amp,2); % 1- filtered raw signal; >1-filtered shuffle signal

if nSurrogates ==1  % MVL for filtered raw signal with one time series
    % Create composite signal
    compositeSig = amp(:).*exp(1i * phase(:));    
    % Mean length of composite signal.
    rawMVL = abs(mean(compositeSig(:))); 
    shuffledMVL = 0;
else % MVL for filtered shuffled signal with nSurrogates
    shuffledMVL = zeros(1, nSurrogates);
    parfor nSurr = 1:nSurrogates        
        %Create composite signal
        compositeSig = amp(:, nSurr).*exp(1i * phase(:, nSurr));        
        %Compute the mean length of composite signal for each Surrogate
        shuffledMVL(1, nSurr) = abs(mean(compositeSig));       
    end
    rawMVL = 0;  
end

end
