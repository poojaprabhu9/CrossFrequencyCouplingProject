folderSourceString = fileparts(pwd); % folder in which programs and data for this project are kept. This line works only if 'pwd' returns the programs folder. Otherwise, simply specify this string, e.g. 'E:\Monkey\MATLAB\CrossFrequencyCouplingProject\';
gridType = 'Microelectrode';

%%%%%%%%%%%%%%%%%%%%%%%%%% Choice of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002';
% monkeyName = 'alpaH'; expDate = '050817'; protocolName = 'GRF_002';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';

rfData = load([monkeyName 'MicroelectrodeRFData.mat']); % selecting good electrodes as per RMS values from Dubey and Ray, Sci Rep, 2020
goodElectrodes = rfData.highRMSElectrodes;

numElectrodes = length(goodElectrodes);

folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);
folderNameMP = fullfile('../data',monkeyName,gridType,expDate,protocolName,'mpAnalysis'); % The path must be relative when using Windows (i.e., should not start with c: etc. So this part only works when we are in the programs folder)
makeDirectory(folderNameMP);

tmp = load(fullfile(folderName,'segmentedData','LFP','lfpInfo.mat'),'timeVals');
timeVals = tmp.timeVals;
Fs = round(1/((timeVals(2)-timeVals(1))));

Max_iterations = 500; % number of iterations

for i=1:numElectrodes
    eNum = goodElectrodes(i);

    tmp = load(fullfile(folderName,'segmentedData','LFP',['elec' num2str(eNum) '.mat']),'analogData');
    inputSignal = tmp.analogData;

    tag = sprintf('elec%d/',eNum);

    % Import the data
    X(:,:,1) = inputSignal';
    L = size(X,1);
    signalRange = [1 L]; % full range
    importData(X,folderNameMP,tag,signalRange,Fs);

    % perform Gabor decomposition
    Numb_points = L; % length of the signal
    runGabor(folderNameMP,tag,Numb_points,Max_iterations);
end
