% runSaveData
clear all;
addpath(genpath('E:\Monkey\MATLAB\Repositories\CommonPrograms-master'))
addpath(genpath('E:\Monkey\MATLAB\Repositories\chronux_2_12'));
folderSourceString = 'E:\Monkey\MATLAB\CrossFrequencyCouplingProject\';

%%%%%%%%%%%%%Choice of parameters%%%%%%%%%%%%%

monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002';
% monkeyName = 'alpaH'; expDate = '050817'; protocolName = 'GRF_002';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';
methodVar= {'mi','plv','glm', 'tort','ozkurt'};
modality='LFP'; %'LFP' or 'ECoG'
electrodeDistanceVal='0'; %'0' '400'
filterMethod={'morlet','FIR','MP'};
removeEvokedResponse=1; tapers = [1 1];
tuning='SF4Ori90'; %'SF' 'Ori' 'SF4Ori90' 'SF2Ori90' 'Size'

saveData(monkeyName,expDate,protocolName,folderSourceString,removeEvokedResponse,tapers, modality, electrodeDistanceVal, tuning, filterMethod, methodVar)

