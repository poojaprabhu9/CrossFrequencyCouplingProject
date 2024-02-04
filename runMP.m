clear all; close all;
% Matching pursuit on the electrodes
addpath(genpath('E:\Monkey\MATLAB\Repositories\CommonPrograms-master'))
addpath(genpath('E:\Monkey\MATLAB\Repositories\chronux_2_12'));
addpath(genpath('E:\Monkey\MATLAB\CrossFrequencyCouplingProject\ForToolbox\MP-master'));

elecNum=90;
folderSourceString = 'E:\Monkey\MATLAB\CrossFrequencyCouplingProject\';
% monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; gridType='Microelectrode';
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; gridType='Microelectrode';
monkeyName = 'alpaH'; expDate = '050817'; protocolName = 'GRF_002'; gridType='Microelectrode';

for e=1:elecNum
        tmp = load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(e) '.mat']),'analogData');
        analogData = tmp.analogData;

        inputSignal=analogData(1,:);
        folderName = '../MPalpaHSize/'; 
        tag = sprintf('elec%d/',e);
        
        % Import the data
        X(:,:,1) = inputSignal';
        L=size(X,1);
        signalRange = [1 L]; % full range
        Fs=2000;
        importData(X,folderName,tag,signalRange,Fs);
        
        % perform Gabor decomposition
        Numb_points = L; % length of the signal
        Max_iterations = 100; % number of iterations
        runGabor(folderName,tag,Numb_points, Max_iterations);
        
        trialNum=size(X,2); % plot Trial number
        rSignal=zeros(size(analogData,1),size(analogData,2));
        for i=1:trialNum
        
            gaborInfo{1} = getGaborData(folderName,tag,1);
            
            
            % Reconstruct signal
            wrap=1;
            atomList=[]; % all atoms
            
            if isempty(atomList)
                disp(['Reconstructing trial ' num2str(i) ', all atoms']);
            else
                disp(['Reconstructing trial ' num2str(i) ', atoms ' num2str(atomList(1)) ':' num2str(atomList(end))]);
            end
            
            rSignal(i,:) = reconstructSignalFromAtomsMPP(gaborInfo{1}{trialNum}.gaborData,L,wrap,atomList);
        end
        filetoSave=fullfile(folderName,tag,['elec' num2str(e) '.mat']);
        save(filetoSave,"rSignal");
        clear gaborInfo
end
