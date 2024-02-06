% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms
% Chronux: http://chronux.org/

folderSourceString = fileparts(pwd); % folder in which programs and data for this project are kept. This line works only if 'pwd' returns the programs folder. Otherwise, simply specify this string, e.g. 'E:\Monkey\MATLAB\CrossFrequencyCouplingProject\';

%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002';
% monkeyName = 'alpaH'; expDate = '050817'; protocolName = 'GRF_002';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';

removeEvokedResponse=1; tapers = [1 1];

tuning='SF4Ori90'; %'SF2Ori90' 'SF4Ori90' 'SF' 'Ori' 'Size'
metrics='pac'; %'others' 'pac' 'pval' 
correction= 'none'; % 'fdr' 'bonf' 'none'
pThresh=0.05; % 0.05,  0.01
methodVar= {'mi','plv','glm', 'tort','ozkurt'};
modality='LFP'; %'LFP' or 'ECoG'
electrodeDistanceVal='0'; % '0' or '400'
filterMethod={'morlet','FIR','MP'};



firingThresh = 10;
snrThresh=1.5;
gridType='Microelectrode';
cutoffs=[firingThresh snrThresh 0 0];

%Choose  the Electrodes based on modality
switch modality

    case 'ECoG'

        if strcmp(monkeyName,'alpaH')

            ecogElectrodes=[82 85 86 88 89];

        elseif strcmp(monkeyName,'kesariH')

            ecogElectrodes=[85 86 88 89];
        end
        lfpElectrodes=ecogElectrodes;
        spikeElectrodes=ecogElectrodes;

    case 'LFP'

        
        [spikeElectrodes,electrodesToUse,firingRate,snr,totalSpikes] = getGoodSpikeElectrodes(monkeyName,expDate,protocolName,folderSourceString,cutoffs);

end

displayData(monkeyName,expDate,protocolName,removeEvokedResponse,tapers,spikeElectrodes,analysisPeriodNum,modality, electrodeDistanceVal, tuning, filterMethod, methodVar, metrics, correction, pThresh)
