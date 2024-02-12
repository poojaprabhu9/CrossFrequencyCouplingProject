function [rawnMVL, shufflednMVL] = methodnMVL (amp, phase)
% normalised mean vector length, Ozkurt 2010

% number of Surrogates
nSurrogates = size(amp,2); % 1- filtered raw signal; >1-filtered shuffle signal

% Number of time (or data) points
nTime = size(amp,1); 

if nSurrogates ==1  % nMVL for filtered raw signal with one time series
    % Create composite signal
    compositeSig = amp(:) .* exp(1i * phase(:));
    % Normalise mean vector length
    rawnMVL = (1./sqrt(nTime)) * abs(mean(compositeSig)) / sqrt(mean(amp(:) .* amp(:))); 
    shufflednMVL = 0;
else % nMVL for filtered shuffled signal with nSurrogates
    shufflednMVL = zeros(1, nSurrogates);
    parfor nSurr = 1:nSurrogates
        % Create composite signal
        compositeSig = amp(:, nSurr).*exp(1i * phase(:, nSurr));
        % Normalise
        shufflednMVL(1, nSurr) = (1./sqrt(nTime)) * abs(mean(compositeSig)) / sqrt(mean(amp(:, nSurr) .* amp(:, nSurr)));
    end
    rawnMVL = 0;  
end

end
