% This is the main program for saving spike-triggered average (STA),
% spike-field coherence (SFC) and field-field coherence (FFC) data. For
% consistency FFC is calculated on the same set of spike-LFP electrodes as
% SFC and STA.

function saveData(monkeyName,expDate,protocolName,folderSourceString,removeEvokedResponse,tapers)

rfData = load([monkeyName 'MicroelectrodeRFData.mat']);

goodElectrodes = rfData.highRMSElectrodes;
lfpElectrodes = goodElectrodes(goodElectrodes<=81);
% goodSpikeElectrodes = getGoodSpikeElectrodes; % Write a function to get
% good spike electrodes. For now, taking all highRMSElectrodes
spikeElectrodes = lfpElectrodes;

saveDataAllElectrodePairs(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,removeEvokedResponse,tapers);

end

function saveDataAllElectrodePairs(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,removeEvokedResponse,tapers)

folderSave = 'savedData';
makeDirectory(folderSave);

analysisPeriodList{1} = [-0.5 0];
analysisPeriodList{2} = [0.25 0.75];

%%%%%%%%%%%%%%%%%%%%%%% Coherence across all pairs %%%%%%%%%%%%%%%%%%%%%%%%
numElectrodes = length(spikeElectrodes);

for i=1:numElectrodes

    spikeElectrode = spikeElectrodes(i);

    disp([num2str(i) ' of ' num2str(numElectrodes) ': ' monkeyName expDate protocolName 'elec' num2str(spikeElectrode)]);

    fileToSave = fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(spikeElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '.mat']);

    if ~exist(fileToSave,'file')

        clear('ffc','ffphi','sfc','sfphi','ffcFreq','sfcFreq','staVals','xsSTA');
        [ffc,ffphi,sfc,sfphi,ffcFreq,sfcFreq,staVals,xsSTA] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,lfpElectrodes,analysisPeriodList,removeEvokedResponse,tapers);

        % save data
        save(fileToSave,'ffc','ffphi','sfc','sfphi','ffcFreq','sfcFreq','staVals','xsSTA');
    else
        disp('file exists');
    end
end
end

function [ffc,ffphi,sfc,sfphi,ffcFreq,sfcFreq,staVals,xsSTA] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,lfpElectrodes,analysisPeriodList,removeEvokedResponse,tapers)

gridType='Microelectrode';
% Get parameter combinations
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat'),'parameterCombinations');
parameterCombinations = tmp.parameterCombinations;

% Get bad trials
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'),'badTrials');
badTrials = tmp.badTrials;

% LFP timeVals
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat'),'timeVals');
timeVals = tmp.timeVals;

numSFs = size(parameterCombinations,4); % including the "all condition"
numOri = size(parameterCombinations,5); % including the "all condition"
numPeriods = length(analysisPeriodList);

Fs = 2000;
for i=1:numPeriods
    xPos{i} = find(timeVals>analysisPeriodList{i}(1),1) + (1:round(Fs*diff(analysisPeriodList{i})));
end

% Get spike and LFP data from the spike electrode
% Get Spike Data
clear spikeData
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','Spikes',['elec' num2str(spikeElectrode) '_SID0.mat']));
spikeData = tmp.spikeData;

% Get LFP data from the spike electrode
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(spikeElectrode) '.mat']),'analogData');
analogDataSpike=tmp.analogData;

numLFPElectrodes = length(lfpElectrodes);

% Set up multitaper
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 1;

for i=1:numLFPElectrodes
    tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(lfpElectrodes(i)) '.mat']),'analogData');
    analogData = tmp.analogData;

    for f=1:numSFs
        for o=1:numOri
            goodPos = setdiff(parameterCombinations{1,1,1,f,o},badTrials);

            for p=1:numPeriods

                spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/Fs);
                lfpSpike = analogDataSpike(goodPos,xPos{p});
                lfp = analogData(goodPos,xPos{p});

                if removeEvokedResponse
                    lfpSpike =removeMeanResponse(lfpSpike);
                    lfp=removeMeanResponse(lfp);
                end

                % FFC and SFC
                [ffc(i,f,o,p,:),ffphi(i,f,o,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
                [sfc(i,f,o,p,:),sfphi(i,f,o,p,:),~,~,~,sfcFreq]=coherencycpb(lfp',spkb,params);

                % Spike-triggered average
                [staVals(i,f,o,p,:),~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);

                % phase-amplitude coupling

            end
        end
    end
end
end

function y=removeMeanResponse(analogData)
y = analogData-repmat(mean(analogData),size(analogData,1),1);
end