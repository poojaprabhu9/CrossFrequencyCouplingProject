function  [rawGLM, shuffledGLM] = methodGLM (amp, phase)
% General linear Model PAC method Penny et al., 2008

% number of Surrogates
nSurrogates = size(amp,2); % 1- filtered raw signal; >1-filtered shuffle signal

if nSurrogates ==1  % GLM for filtered raw signal with one time series
    % Building Matrix of regressors. Note : glmfit adds a column of 1s
    X = [cos(phase(:)), sin(phase(:)) ones(size(phase(:)))];     
    % Fit the GLM    
    [~, ~, stats] = glmfit(X, amp(:), 'normal', 'constant', 'off');     
    % Calculate the variance
    rawGLM = 1- sum(stats.resid.^2)/sum((amp(:) - mean(amp(:))).^2); 
    shuffledGLM = 0;
else % GLM for filtered shuffled signal with nSurrogates
    shuffledGLM = zeros(1, nSurrogates);
    parfor nSurr = 1:nSurrogates        
        % Building Matrix of regressors. Note : glmfit adds a column of 1s
        X = [cos(phase(:, nSurr)), sin(phase(:, nSurr)) ones(size(phase(:, nSurr)))];
        % Fit the GLM    
        [~, ~, stats] = glmfit(X, amp(:, nSurr), 'normal', 'constant', 'off');         
        % Calculate the variance
        shuffledGLM(1, nSurr) = 1- sum(stats.resid.^2)/sum((amp(:, nSurr) - mean(amp(:, nSurr))).^2);       
    end
    rawGLM = 0;  
end

end

