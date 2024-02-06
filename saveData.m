% Function for saving Field-Field Coherence, Spike Field Coherence and
% Spike Triggered Average and PAC using different methods

function saveData(monkeyName,expDate,protocolName,folderSourceString,removeEvokedResponse,tapers,modality,electrodeDistanceVal,sVarName,sPos,oPos,methodVar,filterMethod,nSurr,useMPFlag)

% get good spike eletrodes
firingThresh = 5;
snrThresh = 1.5;

cutoffs = [firingThresh snrThresh 0 0];
disp(['Getting spike electrodes with FR>' num2str(firingThresh) ' and snr>' num2str(snrThresh)]);
spikeElectrodes = getGoodSpikeElectrodes(monkeyName,expDate,protocolName,folderSourceString,cutoffs);
disp([num2str(length(spikeElectrodes)) ' spike electrodes obtained']);

% load electrodes having high RMS value
rfData = load([monkeyName 'MicroelectrodeRFData.mat']); % selecting good electrodes as per RMS values from Dubey and Ray, Sci Rep, 2020
goodElectrodes = rfData.highRMSElectrodes;

% Select good LFP or ECoG electrodes
switch modality
    case 'ECoG'
        if strcmp(monkeyName,'alpaH')
            ecogElectrodes = [82 85 86 88 89];
        elseif strcmp(monkeyName,'kesariH')
            ecogElectrodes = [85 86 88 89];
        end
        lfpElectrodes = ecogElectrodes;
        spikeElectrodes = ecogElectrodes;

    case 'LFP'
        lfpElectrodes = goodElectrodes(goodElectrodes<=81);
end

saveDataEachElectrode(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,removeEvokedResponse,tapers,modality,electrodeDistanceVal,sVarName,sPos,oPos,methodVar,filterMethod,nSurr,useMPFlag);
end

function saveDataEachElectrode(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,removeEvokedResponse,tapers,modality,electrodeDistanceVal,sVarName,sPos,oPos,methodVar,filterMethod,nSurr,useMPFlag)

% Selection of LFP electrodes based on choice of electrode distance between
% LFP and paricular spike electrode
if strcmp(modality,'ECoG')
    electrodeDistanceList{1} = 0;
    electrodeDistanceVal='0';

elseif strcmp(modality,'LFP')
    switch electrodeDistanceVal
        case '0'
            electrodeDistanceList{1} = 0;
        case '400'
            electrodeDistanceList{1} = 400;
        case '400to1600'
            electrodeDistanceList{1} = [400 1600];
        case '1600to2400'
            electrodeDistanceList{1} = [1600 2400];
        case '2400to4000'
            electrodeDistanceList{1} = [2400 4000];
    end
end

folderSave = 'savedData';
makeDirectory(folderSave);

analysisPeriodList{1} = [-0.5 0];
analysisPeriodList{2} = [0.25 0.75];

%%%%%%%%%%%%%%%% Coherence and Coupling across all pairs %%%%%%%%%%%%%%%%%%

% File to save
fileToSave = fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponse) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sPos) '_o' num2str(oPos) '_' methodVar '_' filterMethod '_nSur' num2str(nSurr) '_MP' num2str(useMPFlag) '.mat']);

if exist(fileToSave,'file')
    disp([fileToSave ' exists']);
else
    numElectrodes = length(spikeElectrodes);

    ffc=[]; ffphi=[]; sfc=[]; sfphi=[]; staVals=[]; pacmat=[]; pval=[]; lfpElectrodesToUse=[];
    for i=1:numElectrodes
        spikeElectrode = spikeElectrodes(i);
        [electrodesToCombine,~] = getElectrodeDistance(lfpElectrodes,spikeElectrode,electrodeDistanceList);
        lfpElectrodesSet=lfpElectrodes(electrodesToCombine{1});

        disp([num2str(i) ' of ' num2str(numElectrodes) ': ' monkeyName expDate protocolName 'elec' num2str(spikeElectrode)]);
        [ffcTMP,ffphiTMP,sfcTMP,sfphiTMP,ffcFreq,sfcFreq,staValsTMP,xsSTA,pacmatTMP,pvalTMP,freqVecPhase,freqVecAmp] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,lfpElectrodesSet,analysisPeriodList,removeEvokedResponse,tapers,sVarName,sPos,oPos,methodVar,filterMethod,nSurr,useMPFlag);
        
        % Concatenate
        ffc = cat(1,ffc,ffcTMP);
        ffphi = cat(1,ffphi,ffphiTMP);
        sfc = cat(1,sfc,sfcTMP);
        sfphi = cat(1,sfphi,sfphiTMP);
        staVals = cat(1,staVals,staValsTMP);
        pacmat = cat(1,pacmat,pacmatTMP);
        pval = cat(1,pval,pvalTMP);
        lfpElectrodesToUse = cat(2,lfpElectrodesToUse,lfpElectrodesSet);
    end
    % save data
    lfpElectrodes = lfpElectrodesToUse;
    save(fileToSave,'ffc','ffphi','sfc','sfphi','ffcFreq','sfcFreq','staVals','xsSTA','pacmat','pval','freqVecPhase','freqVecAmp','spikeElectrodes','lfpElectrodes');
end
end

function [electrodesToCombine,electrodeDistances] = getElectrodeDistance(lfpElectrodes,spkElectrode,electrodeDistanceList)

gridType='Microelectrode';
[spkRow,spkColumn,electrodeArray] = electrodePositionOnGrid(spkElectrode,gridType);

electrodeDistances = zeros(1,length(lfpElectrodes));
for i=1:length(lfpElectrodes)
    [lfpRow,lfpColumn] = find(electrodeArray==lfpElectrodes(i));
    electrodeDistances(i) = 400*sqrt((lfpRow-spkRow).^2 + (lfpColumn-spkColumn).^2);
end

%%%%%%%%%%%%%%%%%%%%%%%% Electrodes to combine %%%%%%%%%%%%%%%%%%%%%%%%%%%%
numCombinations = length(electrodeDistanceList);
electrodesToCombine = cell(1,numCombinations);

for i=1:numCombinations
    dRange = electrodeDistanceList{i};
    if length(dRange)==1
        electrodesToCombine{i} = find(electrodeDistances==dRange);
    else
        electrodesToCombine{i} = intersect(find(electrodeDistances>dRange(1)),find(electrodeDistances<=dRange(2)));
    end
end
end

function [ffc,ffphi,sfc,sfphi,ffcFreq,sfcFreq,staVals,xsSTA,pacmat,pval,freqVecPhase,freqVecAmp] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,lfpElectrodes,analysisPeriodList,removeEvokedResponse,tapers,sVarName,sPos,oPos,methodVar,filterMethod,nSurr,useMPFlag)

% PAC analysis
phaseFreqVec = 2:4:80; % Phase frequencies for PAC
ampFreqVec = 2:10:200; % amplitude frequencies for PAC
width = 3; % width of morlet wavelet
epochFrames = 0;
phFiltOrder = 300;
ampFiltOrder = 300; % values taken from Kramer and Eden, JNM, 2013.

gridType='Microelectrode';
% Get parameter combinations
parameterCombinations = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat'));

% Get bad trials
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'),'badTrials');
badTrials = tmp.badTrials;

% Get stimulus positions
if strcmp(sVarName,'sf') % SFOri protocol
    goodPos = setdiff(parameterCombinations.parameterCombinations{1,1,1,sPos,oPos},badTrials);
    if oPos>length(parameterCombinations.oValsUnique)
        disp(['Using SF ' num2str(parameterCombinations.fValsUnique(sPos)) ', OriAll']);
    else
        disp(['Using SF ' num2str(parameterCombinations.fValsUnique(sPos)) ', Ori' num2str(parameterCombinations.oValsUnique(oPos))]);
    end
else
end

% LFP timeVals
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat'),'timeVals');
timeVals = tmp.timeVals;
Fs = round(1/((timeVals(2)-timeVals(1))));

numPeriods = length(analysisPeriodList);
for i=1:numPeriods
    xPos{i} = find(timeVals>analysisPeriodList{i}(1),1) + (1:round(Fs*diff(analysisPeriodList{i})));
end

%%%%%%%%%%%%%% Get spike and LFP data from the spike electrode %%%%%%%%%%%%
% Get Spike Data
clear spikeData
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','Spikes',['elec' num2str(spikeElectrode) '_SID0.mat']));
spikeData = tmp.spikeData;

if useMPFlag
else
    % Get LFP data from the spike electrode
    tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(spikeElectrode) '.mat']),'analogData');
    analogDataSpike=tmp.analogData;
end

numLFPElectrodes = length(lfpElectrodes);

% Set up multitaper
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 250];
params.trialave = 1;

for i=1:numLFPElectrodes

    if useMPFlag
    else
        tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(lfpElectrodes(i)) '.mat']),'analogData');
        analogData = tmp.analogData;
    end

    for p=1:numPeriods

        spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/Fs);
        lfpSpike = analogDataSpike(goodPos,xPos{p});
        lfp = analogData(goodPos,xPos{p});

        if removeEvokedResponse
            lfpSpike = removeMeanResponse(lfpSpike);
            lfp = removeMeanResponse(lfp);
        end

        % FFC and SFC
        [ffc(i,p,:),ffphi(i,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
        [sfc(i,p,:),sfphi(i,p,:),~,~,~,sfcFreq]=coherencycpb(lfp',spkb,params);

        % Spike-triggered average
        [staTMP,~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);
        staVals(i,p,:) = staTMP{1};

        % phase-amplitude coupling
        [pacmat(i,p,:,:), pval(i,p,:,:),freqVecPhase,freqVecAmp] = getPAC(lfpSpike',Fs,filterMethod,methodVar, ...
                lfp',phaseFreqVec,ampFreqVec,width,nSurr,epochFrames,phFiltOrder,ampFiltOrder);
    end
end
end

function y=removeMeanResponse(analogData)
y = analogData-repmat(mean(analogData),size(analogData,1),1);
end