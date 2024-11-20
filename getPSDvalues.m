% Generating the Power Spectral Density (PSD) for each good spike Electrodes for each Spatial Frequency (SF) and Orientation (Ori)
% The generated PSD was used in Figure 1 and Figure 4

clear
folderSourceString = fileparts(pwd);

%%%%%%%%%%%%%%%%%%%%%%%%%% Choice of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';

% SF values and Ori values to compute PSD
% Choice of parameters
% For SFOri protocols, 5 SFs (0.5, 1, 2, 4 and 8 CPD) and 8 orientations
% 0:22.5:157.5) were used. We give the positions of these parameters as
% inputs. For example, sPos = 2, oPos = 5 for sf=2 CPD and ori = 90
% degrees. To select all orientations, choose oPos = 9
% These are the index values for each SF's and Ori
sfPos = [1 2 3 4 5 6]; 
oriPos = [1 2 3 4 5 6 7 8 9];

% Electrode choices
modality = 'LFP'; %'LFP' or 'ECoG'

% Baseline and Stim period
analysisPeriodList{1} = [-0.5 0];
analysisPeriodList{2} = [0.25 0.75];

% Select good LFP or ECoG electrodes
switch modality
    case 'ECoG'
        if strcmp(monkeyName,'alpaH')
            ecogElectrodes = [82 85 86 88 89];
        elseif strcmp(monkeyName,'kesariH')
            ecogElectrodes = [85 86 88 89];
        end      
    spikeElectrodes = ecogElectrodes;
    case 'LFP'
        % get good spike eletrodes 
        firingThresh = 2;
        snrThresh = 1.5;
        
        cutoffs = [firingThresh snrThresh 0 0];
        disp(['Getting spike electrodes with FR>' num2str(firingThresh) ' and snr>' num2str(snrThresh)]);
        [spikeElectrodes,spikeUnitsElectrodes,spikeUnitsSourceID,~,~,~,~] = getGoodSpikeElectrodes(monkeyName,expDate,protocolName,folderSourceString,cutoffs);
        spikeElectrodes = spikeElectrodes(spikeElectrodes<=81);
end


%%% Params
tapers = [2 3]; %As per Dubey and Ray, 2020
movingwin = [0.25 0.025] ;

% PSD for each good Spike electrodes
numElectrodes = length(spikeElectrodes);
for nElec = 1:numElectrodes
    for sPos = 1:length(sfPos)
        for oPos = 1:length(oriPos)
            spikeElectrode = spikeElectrodes(nElec);
        
            gridType='Microelectrode';
            % Get parameter combinations
            parameterCombinations = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat'));
            
            % Get bad trials
            tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'),'badTrials');
            badTrials = tmp.badTrials;        
            goodPos = setdiff(parameterCombinations.parameterCombinations{1,1,1,sfPos(sPos),oriPos(oPos)},badTrials);
        
            % LFP timeVals
            tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat'),'timeVals');
            timeVals = tmp.timeVals;
            samplingFreq = round(1/((timeVals(2)-timeVals(1))));
            
            numPeriods = length(analysisPeriodList);
            for i=1:numPeriods
                xPos{i} = find(timeVals>analysisPeriodList{i}(1),1) + (1:round(samplingFreq*diff(analysisPeriodList{i})));
            end
        
            tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(spikeElectrode) '.mat']),'analogData');
            analogData=tmp.analogData;
        
            % Set up multitaper
            params.tapers   = tapers;
            params.pad      = -1;
            params.Fs       = samplingFreq;
            params.fpass    = [0 250];
            params.trialave = 1;
        
            for p=1:numPeriods    
                tmpSignal = analogData(goodPos,xPos{p})';
                [energyValues(nElec,p,sPos,oPos,:),mtFreq]=mtspectrumc(tmpSignal,params);
            end
        end
    end
end
% final outcome saved in savedData