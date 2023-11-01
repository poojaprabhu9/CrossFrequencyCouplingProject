
function displayCoherenceData(monkeyName,expDate,protocolName,removeEvokedResponse,tapers,spikeElectrode,f,o,analysisPeriodNum)

rfData = load([monkeyName 'MicroelectrodeRFData.mat']);

goodElectrodes = rfData.highRMSElectrodes;
lfpElectrodes = goodElectrodes(goodElectrodes<=81);
folderSave = 'coherenceData';

electrodeDistanceList{1} = 0;
electrodeDistanceList{2} = 400;
electrodeDistanceList{3} = [400 1600];
electrodeDistanceList{4} = [1600 2400];
electrodeDistanceList{5} = [2400 4000];

distanceColorNames = ['k';'r';'b';'c';'g'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
tmp = load(fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(spikeElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '.mat']));
ffc = squeeze(tmp.ffc(:,f,o,analysisPeriodNum,:));
sfc = squeeze(tmp.sfc(:,f,o,analysisPeriodNum,:));
ffphi = squeeze(tmp.ffphi(:,f,o,analysisPeriodNum,:));
sfphi = squeeze(tmp.sfphi(:,f,o,analysisPeriodNum,:));
sta = cell2mat(squeeze(tmp.staVals(:,f,o,analysisPeriodNum)));

ffcFreq = tmp.ffcFreq;
sfcFreq = tmp.sfcFreq;
xsSTA = tmp.xsSTA;

electrodeDistances = getElectrodeDistance(lfpElectrodes,spikeElectrode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numFreqs = length(ffcFreq);
numCombinations = length(electrodeDistanceList);
meanFFC = zeros(numCombinations,numFreqs);
meanSFC = zeros(numCombinations,numFreqs);
meanFFPhi = zeros(numCombinations,numFreqs);
meanSFPhi = zeros(numCombinations,numFreqs);
meanSTA = zeros(numCombinations,length(xsSTA));

numElectrodes = zeros(1,numCombinations);

for i=1:numCombinations
    dRange = electrodeDistanceList{i};
    if length(dRange)==1
        electrodesToCombine = find(electrodeDistances==dRange);
    else
        electrodesToCombine = intersect(find(electrodeDistances>dRange(1)),find(electrodeDistances<=dRange(2)));
    end
    
    numElectrodes(i) = length(electrodesToCombine);
    if numElectrodes(i)>0
        meanFFC(i,:) = mean(ffc(electrodesToCombine,:),1);
        meanSFC(i,:) = mean(sfc(electrodesToCombine,:),1);
        
        meanFFPhi(i,:) = circ_mean(ffphi(electrodesToCombine,:),[],1);
        meanSFPhi(i,:) = circ_mean(sfphi(electrodesToCombine,:),[],1);
        
        meanSTA(i,:) = mean(sta(electrodesToCombine,:),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:numCombinations
    subplot(221)
    plot(ffcFreq,meanFFC(i,:),distanceColorNames(i));
    hold on;

    subplot(222)
    plot(sfcFreq,meanSFC(i,:),distanceColorNames(i));
    hold on;

    subplot(223)
    plot(xsSTA,meanSTA(i,:),distanceColorNames(i));
    hold on;
end

end
function distance = getElectrodeDistance(lfpElectrodes,spkElectrode)
gridType='Microelectrode';
[spkRow,spkColumn,electrodeArray] = electrodePositionOnGrid(spkElectrode,gridType);

distance = zeros(1,length(lfpElectrodes));
for i=1:length(lfpElectrodes)
    [lfpRow,lfpColumn] = find(electrodeArray==lfpElectrodes(i));
    distance(i) = 400*sqrt((lfpRow-spkRow).^2 + (lfpColumn-spkColumn).^2);
end
end