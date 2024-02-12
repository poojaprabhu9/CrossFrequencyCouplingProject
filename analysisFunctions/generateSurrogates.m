function [shuffledAmp, shuffledPhase] = generateSurrogates (rawAmp, rawPhase, nSurrogates) 

% Shuffle time indices to generate shuffled amplitude and phase signals
% Set seed 
rng('default');
rng(1);

% number of time points
nTime = size(rawAmp{1,1},1);

% randomise indices to generate shuffled amplitude signals
randomIdx = randi([1 nTime],nSurrogates,nTime);
disp(['Generating ' num2str(nSurrogates) ' shuffled amplitude signals']);

nAmpBin = size(rawAmp,2);
shuffledAmp = cell(1,nAmpBin);

for ampBin = 1:nAmpBin
    for nSurr = 1:nSurrogates
        shuffledAmp{1,ampBin}(:,nSurr) = rawAmp{1,ampBin}(randomIdx(nSurr,:),:);
    end
end


% Set seed 
rng('default');
rng(2);    

% Randomise indices to generate shuffled phase signals
randomIdx= randi([1 nTime],nSurrogates,nTime);
disp(['Generating ' num2str(nSurrogates) ' shuffled phase signals']);

nPhaseBin=size(rawPhase,2);
shuffledPhase = cell(1,nPhaseBin); 

for phaseBin=1:nPhaseBin
    for nSurr=1:nSurrogates
        shuffledPhase{1,phaseBin}(:,nSurr) = rawPhase{1,phaseBin}(randomIdx(nSurr,:),:);
    end
end

end