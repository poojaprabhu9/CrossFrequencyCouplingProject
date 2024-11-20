% Function for saving Field-Field Coherence, Spike Field Coherence and
% Spike Triggered Average and PAC using different methods

function saveData(monkeyName,expDate,protocolName,folderSourceString,removeEvokedResponseFlag,tapers,modality,electrodeDistanceVal,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag)

% get good spike eletrodes
firingThresh = 2;
snrThresh = 1.5;

cutoffs = [firingThresh snrThresh 0 0];
disp(['Getting spike electrodes with FR>' num2str(firingThresh) ' and snr>' num2str(snrThresh)]);
[spikeElectrodes,spikeUnitsElectrodes,spikeUnitsSourceID,~,~,~,~] = getGoodSpikeElectrodes(monkeyName,expDate,protocolName,folderSourceString,cutoffs);
spikeElectrodes = spikeElectrodes(spikeElectrodes<=81);
spikeUnitsSourceID(spikeUnitsElectrodes>81)= [];
spikeUnitsElectrodes(spikeUnitsElectrodes>81) = [];
disp([num2str(length(spikeElectrodes)) ' spike electrodes obtained']);
disp([num2str(length(spikeUnitsElectrodes)) ' spiking units obtained']);

% getting corresponding SID number for each spike electrode
SIDnum = cell(1,length(spikeElectrodes));
for nSpikeElec = 1:length(spikeElectrodes)
    SIDnum{nSpikeElec} = spikeUnitsSourceID(spikeUnitsElectrodes==spikeElectrodes(nSpikeElec));
end

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
        SIDnum = cell(1,length(spikeElectrodes));
        for nSpikeElec = 1:length(spikeElectrodes)
            SIDnum{nSpikeElec} = 1; % just for the pipeline, however we do not consider spike data for ECoG
        end        

    case 'LFP'
        lfpElectrodes = goodElectrodes(goodElectrodes<=81);
end

saveDataEachElectrode(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,SIDnum,removeEvokedResponseFlag,tapers,modality,electrodeDistanceVal,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag);
end

function saveDataEachElectrode(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,SIDnum,removeEvokedResponseFlag,tapers,modality,electrodeDistanceVal,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag)

% Selection of LFP electrodes based on choice of electrode distance between
% LFP and particular spike electrode
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
fileToSave = fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sPos) '_o' num2str(oPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']);

if exist(fileToSave,'file')
    disp([fileToSave ' exists']);
else
    numElectrodes = length(spikeElectrodes);

    ffc=[]; ffphi=[]; sfc=[]; sfphi=[]; staVals=[]; pac=[]; anglePac=[]; normalisedPac=[]; meanAmp=[]; surrogatePac=[]; tval=[]; pval=[]; lfpElectrodesToUse=[];
    for i=1:numElectrodes
        spikeElectrode = spikeElectrodes(i);
        SID = SIDnum{i};
        [electrodesToCombine,~] = getElectrodeDistance(lfpElectrodes,spikeElectrode,electrodeDistanceList);
        lfpElectrodesSet=lfpElectrodes(electrodesToCombine{1});

        disp([num2str(i) ' of ' num2str(numElectrodes) ': ' monkeyName expDate protocolName 'elec' num2str(spikeElectrode)]);
        [ffcTMP,ffphiTMP,sfcTMP,sfphiTMP,ffcFreq,sfcFreq,staValsTMP,xsSTA,pacTMP,anglePacTMP,normalisedPacTMP,meanAmpTMP,surrogatePacTMP,tvalTMP,pvalTMP,centerAmpFreq,centerPhaseFreq] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,SID,lfpElectrodesSet,analysisPeriodList,removeEvokedResponseFlag,tapers,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag);
        
        % Concatenate
        ffc = cat(1,ffc,ffcTMP);
        ffphi = cat(1,ffphi,ffphiTMP);
        sfc = cat(1,sfc,sfcTMP);
        sfphi = cat(1,sfphi,sfphiTMP);
        staVals = cat(1,staVals,staValsTMP);
        pac = cat(1,pac,pacTMP);
        anglePac = cat(1,anglePac,anglePacTMP);
        normalisedPac = cat(1,normalisedPac,normalisedPacTMP);
        surrogatePac = cat(1,surrogatePac,surrogatePacTMP);
        meanAmp = cat(1,meanAmp,meanAmpTMP);
        tval = cat(1,tval,tvalTMP);
        pval = cat(1,pval,pvalTMP);
        lfpElectrodesToUse = cat(2,lfpElectrodesToUse,lfpElectrodesSet);
    end
    % save data
    lfpElectrodes = lfpElectrodesToUse;
    save(fileToSave,'ffc','ffphi','sfc','sfphi','ffcFreq','sfcFreq','staVals','xsSTA','pac','anglePac','normalisedPac','meanAmp','surrogatePac','tval','pval','centerAmpFreq','centerPhaseFreq','spikeElectrodes','lfpElectrodes');
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

function [ffc,ffphi,sfc,sfphi,ffcFreq,sfcFreq,staVals,xsSTA,pac,anglePac,normalisedPac,meanAmp,surrogatePac,tval,pval,centerAmpFreq,centerPhaseFreq] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,SID,lfpElectrodes,analysisPeriodList,removeEvokedResponseFlag,tapers,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag)

% Parameters for PAC analysis
phaseFreq = 2:4:70; % Phase frequencies for PAC
ampFreq = 2:10:500; % amplitude frequencies for PAC
alpha = 0.05; % confidence interval for ttest (Surrogate v/s raw PAC significance) 

% Get parameter combinations
gridType='Microelectrode';
parameterCombinations = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat'));

% Get bad trials
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'),'badTrials');
badTrials = tmp.badTrials;

% Get stimulus positions
if strcmp(sVarName,'sf') % SFOri protocol
    goodPos = setdiff(parameterCombinations.parameterCombinations{1,1,1,sPos,oPos},badTrials);
    if (oPos>length(parameterCombinations.oValsUnique)) && (sPos>length(parameterCombinations.fValsUnique))
        disp(['Using SFAll and OriAll']);
    elseif oPos>length(parameterCombinations.oValsUnique)
        disp(['Using SF ' num2str(parameterCombinations.fValsUnique(sPos)) ', OriAll']);
    else
        if sPos>length(parameterCombinations.fValsUnique)
            disp(['Using SFAll, Ori ' num2str(parameterCombinations.oValsUnique(oPos))]);
        else
            disp(['Using SF ' num2str(parameterCombinations.fValsUnique(sPos)) ', Ori' num2str(parameterCombinations.oValsUnique(oPos))]);
        end
    end
elseif strcmp(sVarName,'sz') % SizeOri protocol
    goodPos = setdiff(parameterCombinations.parameterCombinations{1,1,sPos,1,oPos},badTrials);
    if oPos>length(parameterCombinations.oValsUnique)
        disp(['Using Size ' num2str(parameterCombinations.sValsUnique(sPos)*3) ', OriAll']);
    else
        if sPos>length(parameterCombinations.sValsUnique)
            disp(['Using SizeAll, Ori ' num2str(parameterCombinations.oValsUnique(oPos))]);
        else
            disp(['Using Size ' num2str(parameterCombinations.sValsUnique(sPos)*3) ', Ori' num2str(parameterCombinations.oValsUnique(oPos))]);
        end
    end
end

% LFP timeVals
tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat'),'timeVals');
timeVals = tmp.timeVals;
samplingFreq = round(1/((timeVals(2)-timeVals(1))));

numPeriods = length(analysisPeriodList);
for i=1:numPeriods
    xPos{i} = find(timeVals>analysisPeriodList{i}(1),1) + (1:round(samplingFreq*diff(analysisPeriodList{i})));
end

%%%%%%%%%%%%%% LFP data from the spike electrode %%%%%%%%%%%%
if useMPFlag
    % Get LFP data from the spike electrode
    tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','lfpMP',['elec' num2str(spikeElectrode) '.mat']),'analogData');
    analogDataSpike=tmp.analogData;
else
    % Get LFP data from the spike electrode
    tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(spikeElectrode) '.mat']),'analogData');
    analogDataSpike=tmp.analogData;
    analogDataSpike = removeMean(analogDataSpike);
end

numLFPElectrodes = length(lfpElectrodes);

% Set up multitaper
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = samplingFreq;
params.fpass    = [0 250];
params.trialave = 1;

for i=1:numLFPElectrodes

    if useMPFlag
        tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','lfpMP',['elec' num2str(lfpElectrodes(i)) '.mat']),'analogData');
        analogData = tmp.analogData;
    else
        tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(lfpElectrodes(i)) '.mat']),'analogData');
        analogData = tmp.analogData;
        analogData = removeMean(analogData);
    end

    for p=1:numPeriods

        % Get SFC and STA for each spiking unit corresponding to each spike electrode
        clear sfcPerUnit sfphiPerUnit
        for s = 1:length(SID)
            % Get Spike Data
            clear spikeData
            tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','sortedSpikes',['elec' num2str(spikeElectrode) '_SID' num2str(SID(s)) '.mat']));
            spikeData = tmp.spikeData;
    
            spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/samplingFreq);
            lfpSpike = analogDataSpike(goodPos,xPos{p});
            lfp = analogData(goodPos,xPos{p});
            
            % Use the entire time series for PAC to get rid off edge
            % effects after filtering
            lfpSpikeAll = analogDataSpike(goodPos,:);
            lfpAll = analogData(goodPos,:);
    
            if removeEvokedResponseFlag
                lfpSpike = removeEvokedResponse(lfpSpike);
                lfp = removeEvokedResponse(lfp);
                lfpSpikeAll = removeEvokedResponse(lfpSpikeAll);
                lfpAll = removeEvokedResponse(lfpAll);
            end

            % SFC for each unit
            [sfcPerUnit(s,:),sfphiPerUnit(s,:),~,~,~,sfcFreq] = coherencycpb(lfp',spkb,params);
        
             % Spike-triggered average(STA) for each unit
            [staTMP,~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);
            if ~isempty(staTMP{1})
                staValsPerUnit(s,:) = staTMP{1};
            else
                staValsPerUnit(s,:) = zeros(size(xsSTA));
            end
        end

        % FFC and SFC
        [ffc(i,p,:),ffphi(i,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
        sfc(i,p,:) = mean(sfcPerUnit,1);
        sfphi(i,p,:) = mean(sfphiPerUnit,1);

        %STA
        staVals(i,p,:) = mean(staValsPerUnit,1);

        % phase-amplitude coupling 
        [pac(i,p,:,:),anglePac(i,p,:,:),normalisedPac(i,p,:,:),meanAmp(i,p,:,:,:),surrogatePac(i,p,:,:,:),tval(i,p,:,:),pval(i,p,:,:),centerAmpFreq,centerPhaseFreq] = getPAC (lfpSpikeAll',ampFreq,lfpAll',phaseFreq,samplingFreq,filterName, ...
                                    pacMethod,nSurrogates,alpha,xPos{p});
    end
end
end

function y=removeEvokedResponse(analogData)
y = analogData-repmat(mean(analogData),size(analogData,1),1);
end

function y=removeMean(analogData)
y = analogData-repmat(mean(analogData,2),1,size(analogData,2));
end