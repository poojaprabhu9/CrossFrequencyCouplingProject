% computing and saving Spike field coherence (SFC), Field-field coherence (FFC), spike
% Triggered Average (STA) and Phse-Amplitude coupling (PAC)

% Executes for one monkey at a time

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms
% Chronux: http://chronux.org/

clear;clf;
folderSourceString = fileparts(pwd); % folder in which programs and data for this project are kept. This line works only if 'pwd' returns the programs folder. Otherwise, simply specify this string, e.g. 'E:\Monkey\MATLAB\CrossFrequencyCouplingProject\';

%%%%%%%%%%%%%%%%%%%%%%%%%% Choice of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Protocol details
% Comment one of the below line
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';

% Electrode choices
modality = 'LFP'; %'LFP' or 'ECoG'
electrodeDistanceVal = '0'; %'0' '400'

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];

% Choice of parameters
% For SFOri protocols, 5 SFs (0.5, 1, 2, 4 and 8 CPD) and 8 orientations
% 0:22.5:157.5) were used. We give the positions of these parameters as
% inputs. For example, sPos = 2, oPos = 5 for sf=2 CPD and ori = 90
% degrees. To select all orientations, choose oPos = 9;

% Protocol Name: SFOri
sVarName = 'sf';  
sPos = 4; oPos = 5;

% Methods used for PAC analysis
pacMethod = 'klmi'; % 'mvl','plv','glm', 'klmi','nmvl';
filterName = 'fir'; % 'wavelet','fir', 'butter';
nSurrogates = 0; % Number of surrogates
useMPFlag = 0;

saveData(monkeyName,expDate,protocolName,folderSourceString,removeEvokedResponseFlag,tapers,modality,electrodeDistanceVal,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag);
