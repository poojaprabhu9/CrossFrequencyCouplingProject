%function for saving different methods of Phase Amplitude Coupling,
%Field-Field Coherence, Spike Field Coherence and Spike Triggered
%Average
function saveData(monkeyName,expDate,protocolName,folderSourceString,removeEvokedResponse,tapers, modality, electrodeDistanceVal, tuning, filterMethod, methodVar)

% get good spike eletrodes

firingThresh = 10;
snrThresh=1.5;

cutoffs=[firingThresh snrThresh 0 0];

[spikeElectrodes,~,~,~,~] = getGoodSpikeElectrodes(monkeyName,expDate,protocolName,folderSourceString,cutoffs);

%load electrodes having high RMS value
rfData = load([monkeyName 'MicroelectrodeRFData.mat']);
%selecting good electrodes as per RMS values from Dubey and Ray, Sci Rep, 2020 
goodElectrodes = rfData.highRMSElectrodes;

% Select good LFP or ECoG electrodes 
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

        lfpElectrodes = goodElectrodes(goodElectrodes<=81);
end

   
saveDataEachElectrode(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,removeEvokedResponse,tapers, modality, electrodeDistanceVal, tuning, filterMethod, methodVar);

end


function saveDataEachElectrode(monkeyName,expDate,protocolName,folderSourceString,lfpElectrodes,spikeElectrodes,removeEvokedResponse,tapers,modality, electrodeDistanceVal, tuning, filterMethod, methodVar)

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

%%%%%%%%%%%%%%%%%%%%%%% Coherence and Coupling across all pairs %%%%%%%%%%%%%%%%%%%%%%%%
numElectrodes = length(spikeElectrodes);
% Coherence and PAC between a spike electrode and LFP electrode as per
% choice of electrode distance
for i=1:numElectrodes

    spikeElectrode = spikeElectrodes(i);

    [electrodesToCombine,~] = getElectrodeDistance(lfpElectrodes,spikeElectrode,electrodeDistanceList);

    lfpElectrodesSet=lfpElectrodes(electrodesToCombine{1});

    disp([num2str(i) ' of ' num2str(numElectrodes) ': ' monkeyName expDate protocolName 'elec' num2str(spikeElectrode)]);

    fileToSave = fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(spikeElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_' tuning 'tuning_Allmethods_for_distance' electrodeDistanceVal '.mat']);

    if ~exist(fileToSave,'file')


        clear('ffc','ffphi','sfc','sfphi','ffcFreq','sfcFreq','staVals','xsSTA','pacmat','pval', 'freqvec_ph', 'freqvec_amp');
        [ffc,ffphi,sfc,sfphi,ffcFreq,sfcFreq,staVals,xsSTA,pacmat,pval,freqvec_ph, freqvec_amp] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,lfpElectrodesSet,analysisPeriodList,removeEvokedResponse,tapers, tuning, filterMethod, methodVar);
        % save data
        save(fileToSave,'ffc','ffphi','sfc','sfphi','ffcFreq','sfcFreq','staVals','xsSTA','pacmat','pval','freqvec_ph', 'freqvec_amp');

    else
        disp('file exists');
    end
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



function [ffc,ffphi,sfc,sfphi,ffcFreq,sfcFreq,staVals,xsSTA,pacmat,pval,freqvec_ph, freqvec_amp] = getDataSingleElectrode(monkeyName,expDate,protocolName,folderSourceString,spikeElectrode,lfpElectrodes,analysisPeriodList,removeEvokedResponse,tapers, tuning, filterMethod, methodVar)

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

% number of parameter combinations based on tuning parameter
if strcmp(tuning,'SF')
    numSFs = size(parameterCombinations,4)-1; %exclude 'all'
    numOri = 1;
elseif strcmp(tuning,'SF4Ori90') % for alpaH
    numSFs = 1;
    numOri = 1; 
elseif strcmp(tuning,'SF2Ori90') % for kesariH
    numSFs = 1;
    numOri = 1; 
elseif strcmp(tuning,'Ori')
    numSFs = 1;
    numOri = size(parameterCombinations,5)-1; % %exclude 'all'
elseif strcmp(tuning,'Size')
    numSFs = size(parameterCombinations,3)-1; % %exclude 'all'
    numOri = 1;
end

numPeriods = length(analysisPeriodList);


Fs=2000;
% Phase frequencies for PAC
ph_freq_vec=2:4:80; 
% amplitude frequencies for PAC
amp_freq_vec=2:10:200; 
width=3; %width of morlet wavelet
n_surr=200; %number of surrogates
epochframes=0;
phFiltorder=300; 
ampFiltorder=300; %values taken from Kramer and Eden, JNM, 2013. 


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

    if strcmp(tuning,'SF')
        for f=1:numSFs
            o=5; %90 degrees
%             orientation={'0','22','45','67','90','112','135','157','all'};


            goodPos = setdiff(parameterCombinations{1,1,1,f,o},badTrials);
            disp(['LFP_elec' num2str(i)  ':SF ' num2str(f) 'of' 'Ori' num2str(o)]);

            for p=1:numPeriods

                spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/Fs);
                lfpSpike = analogDataSpike(goodPos,xPos{p});
                lfp = analogData(goodPos,xPos{p});


                if removeEvokedResponse
                    lfpSpike =removeMeanResponse(lfpSpike);
                    lfp=removeMeanResponse(lfp);
                end

                % FFC and SFC
                [ffc(i,f,1,p,:),ffphi(i,f,1,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
                [sfc(i,f,1,p,:),sfphi(i,f,1,p,:),~,~,~,sfcFreq]=coherencycpb(lfp',spkb,params);

                % Spike-triggered average
                [staVals(i,f,1,p,:),~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);

                % phase-amplitude coupling
                nFilterMethod=length(filterMethod);
%                 methods={'esc','mi','plv','tort','ozkurt'};
                nMethods=length(methodVar);
                for fi=1:nFilterMethod
                    filterName=filterMethod{fi};
                    if strcmp(filterName,'MP')

                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(spikeElectrode)],['elec' num2str(spikeElectrode) '.mat']));
                        analogDataSpikeMP=tmpMP.rSignal;
                        lfpSpikeMP = analogDataSpikeMP(goodPos,xPos{p});
                        

                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(lfpElectrodes(i))],['elec' num2str(lfpElectrodes(i)) '.mat']));
                        analogDataMP=tmpMP.rSignal;
                        lfpMP = analogDataMP(goodPos,xPos{p});

                        if removeEvokedResponse
                            lfpSpikeMP =removeMeanResponse(lfpSpikeMP);
                            lfpMP=removeMeanResponse(lfpMP);
                        end

                        for m=1:nMethods
                            measure=methodVar{m};                        
                            [pacmat(i,f,1,p,fi,m,:,:), pval(i,f,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpikeMP', Fs, filterName,measure, ...
                                    lfpMP', ph_freq_vec, amp_freq_vec, width, n_surr, epochframes, phFiltorder,ampFiltorder);

                        end

                    else
                        for m=1:nMethods
                            measure=methodVar{m};                        
                            [pacmat(i,f,1,p,fi,m,:,:), pval(i,f,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpike', Fs, filterName,measure, ...
                                    lfp', ph_freq_vec, amp_freq_vec,width, n_surr, epochframes, phFiltorder,ampFiltorder);
    
                        end
                    end
                end

              disp([measure ' for ' filterName 'for p' num2str(p)]);
            end
        end
    
    elseif strcmp(tuning,'Size')
        for s=1:numSFs
            o=5; %90 degrees
%             orientation={'0','22','45','67','90','112','135','157','all'};


            goodPos = setdiff(parameterCombinations{1,1,s,1,o},badTrials);
            disp(['LFP_elec' num2str(i)  ':Size ' num2str(s) 'of' 'Ori' num2str(o)]);

            for p=1:numPeriods

                spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/Fs);
                lfpSpike = analogDataSpike(goodPos,xPos{p});
                lfp = analogData(goodPos,xPos{p});

                if removeEvokedResponse
                    lfpSpike =removeMeanResponse(lfpSpike);
                    lfp=removeMeanResponse(lfp);
                end

                % FFC and SFC
                [ffc(i,s,1,p,:),ffphi(i,s,1,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
                [sfc(i,s,1,p,:),sfphi(i,s,1,p,:),~,~,~,sfcFreq]=coherencycpb(lfp',spkb,params);

                % Spike-triggered average
                [staVals(i,s,1,p,:),~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);

                % phase-amplitude coupling for each pair of electrodes
                % for different filtering methods and PAC measures
                nFilterMethod=length(filterMethod);
                nMethods=length(methodVar);
                for fi=1:nFilterMethod
                    filterName=filterMethod{fi};
                    if strcmp(filterName,'MP')
                        %for filtering method Matching Pursuit, Considered
                        %the data which is filetered behorehand

                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(spikeElectrode)],['elec' num2str(spikeElectrode) '.mat']));
                        analogDataSpikeMP=tmpMP.rSignal;
                        lfpSpikeMP = analogDataSpikeMP(goodPos,xPos{p});
                        

                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(lfpElectrodes(i))],['elec' num2str(lfpElectrodes(i)) '.mat']));
                        analogDataMP=tmpMP.rSignal;
                        lfpMP = analogDataMP(goodPos,xPos{p});

                        if removeEvokedResponse
                            lfpSpikeMP =removeMeanResponse(lfpSpikeMP);
                            lfpMP=removeMeanResponse(lfpMP);
                        end

                        for m=1:nMethods
                            measure=methodVar{m}; 
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,s,1,p,fi,m,:,:), pval(i,s,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpikeMP', Fs, filterName,measure, ...
                                    lfpMP', ph_freq_vec, amp_freq_vec,  width, n_surr, epochframes, phFiltorder,ampFiltorder);

                        end

                    else
                        for m=1:nMethods
                            measure=methodVar{m};  
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,s,1,p,fi,m,:,:), pval(i,s,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpike', Fs, filterName,measure, ...
                                    lfp', ph_freq_vec, amp_freq_vec, width, n_surr, epochframes, phFiltorder,ampFiltorder);
    
                        end
                    end
                end

              disp([measure ' for ' filterName 'for p' num2str(p)]);
            end
        end


    
    elseif strcmp(tuning,'SF4Ori90')
                   
%             orientation={'0','22','45','67','90','112','135','157','all'};
%             spat_freq={'05','1','2','4','8','all'};
            f=numSFs;o=numOri;


            goodPos = setdiff(parameterCombinations{1,1,1,4,5},badTrials);
            disp(['LFP_elec' num2str(i)  ':SF ' num2str(f) 'of' 'Ori' num2str(o)]);

            for p=1:numPeriods

                spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/Fs);
                lfpSpike = analogDataSpike(goodPos,xPos{p});
                lfp = analogData(goodPos,xPos{p});

                if removeEvokedResponse
                    lfpSpike =removeMeanResponse(lfpSpike);
                    lfp=removeMeanResponse(lfp);
                end

                % FFC and SFC
                [ffc(i,1,1,p,:),ffphi(i,1,1,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
                [sfc(i,1,1,p,:),sfphi(i,1,1,p,:),~,~,~,sfcFreq]=coherencycpb(lfp',spkb,params);

                % Spike-triggered average
                [staVals(i,1,1,p,:),~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);

                % phase-amplitude coupling between electrode pair for each
                % filetering and PAC measures
                nFilterMethod=length(filterMethod);
%                 methods={'esc','mi','plv','tort','ozkurt'};
                nMethods=length(methodVar);
                for fi=1:nFilterMethod
                    filterName=filterMethod{fi};
                    if strcmp(filterName,'MP')
                        % used the already filetered data using Matching
                        % pursuit
                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(spikeElectrode)],['elec' num2str(spikeElectrode) '.mat']));
                        analogDataSpikeMP=tmpMP.rSignal;
                        lfpSpikeMP = analogDataSpikeMP(goodPos,xPos{p});
                        

                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(lfpElectrodes(i))],['elec' num2str(lfpElectrodes(i)) '.mat']));
                        analogDataMP=tmpMP.rSignal;
                        lfpMP = analogDataMP(goodPos,xPos{p});

                        if removeEvokedResponse
                            lfpSpikeMP =removeMeanResponse(lfpSpikeMP);
                            lfpMP=removeMeanResponse(lfpMP);
                        end

                        for m=1:nMethods
                            measure=methodVar{m}; 
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,1,1,p,fi,m,:,:), pval(i,1,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpikeMP', Fs, filterName,measure, ...
                                    lfpMP', ph_freq_vec, amp_freq_vec,  width, n_surr, epochframes, phFiltorder,ampFiltorder);

                        end

                    else
                        for m=1:nMethods
                            measure=methodVar{m}; 
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,1,1,p,fi,m,:,:), pval(i,1,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpike', Fs, filterName,measure, ...
                                    lfp', ph_freq_vec, amp_freq_vec,  width, n_surr, epochframes, phFiltorder,ampFiltorder);
    
                        end
                    end
                end
              disp([measure ' for ' filterName 'for p' num2str(p)]);

            end

     elseif strcmp(tuning,'SF2Ori90')
% for kesariH                   
%             orientation={'0','22','45','67','90','112','135','157','all'};
%             spat_freq={'05','1','2','4','8','all'};
            f=numSFs;o=numOri;


            goodPos = setdiff(parameterCombinations{1,1,1,3,5},badTrials);
            disp(['LFP_elec' num2str(i)  ':SF ' num2str(f) 'of' 'Ori' num2str(o)]);

            for p=1:numPeriods

                spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/Fs);
                lfpSpike = analogDataSpike(goodPos,xPos{p});
                lfp = analogData(goodPos,xPos{p});


                if removeEvokedResponse
                    lfpSpike =removeMeanResponse(lfpSpike);
                    lfp=removeMeanResponse(lfp);
                end

                % FFC and SFC
                [ffc(i,1,1,p,:),ffphi(i,1,1,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
                [sfc(i,1,1,p,:),sfphi(i,1,1,p,:),~,~,~,sfcFreq]=coherencycpb(lfp',spkb,params);

                % Spike-triggered average
                [staVals(i,1,1,p,:),~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);

                % phase-amplitude coupling between electrode pairs for each
                % filtering and PAC measures
                nFilterMethod=length(filterMethod);
%                 methods={'esc','mi','plv','tort','ozkurt'};
                nMethods=length(methodVar);
                for fi=1:nFilterMethod
                    filterName=filterMethod{fi};
                    if strcmp(filterName,'MP')
                        % Use the already filetered data using matching
                        % pursuit
                        
                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(spikeElectrode)],['elec' num2str(spikeElectrode) '.mat']));
                        analogDataSpikeMP=tmpMP.rSignal;
                        lfpSpikeMP = analogDataSpikeMP(goodPos,xPos{p});
                        

                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(lfpElectrodes(i))],['elec' num2str(lfpElectrodes(i)) '.mat']));
                        analogDataMP=tmpMP.rSignal;
                        lfpMP = analogDataMP(goodPos,xPos{p});

                        if removeEvokedResponse
                            lfpSpikeMP =removeMeanResponse(lfpSpikeMP);
                            lfpMP=removeMeanResponse(lfpMP);
                        end

                        for m=1:nMethods
                            measure=methodVar{m}; 
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,1,1,p,fi,m,:,:), pval(i,1,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpikeMP', Fs, filterName,measure, ...
                                    lfpMP', ph_freq_vec, amp_freq_vec,  width, n_surr, epochframes, phFiltorder,ampFiltorder);

                        end

                    else
                        for m=1:nMethods
                            measure=methodVar{m}; 
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,1,1,p,fi,m,:,:), pval(i,1,1,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpike', Fs, filterName,measure, ...
                                    lfp', ph_freq_vec, amp_freq_vec,  width, n_surr, epochframes, phFiltorder,ampFiltorder);
    
                        end
                    end
                end
              disp([measure ' for ' filterName 'for p' num2str(p)]);

            end


     elseif strcmp(tuning,'Ori')
         for o=1:numOri
             if strcmp(monkeyName,'alpaH')
                f=4; %4 cpd
%             spat_freq={'05','1','2','4','8','all'};
             elseif strcmp(monkeyName,'kesariH')
                 f=3;
             end
             goodPos = setdiff(parameterCombinations{1,1,1,f,o},badTrials);
            disp(['LFP_elec' num2str(i)  ':SF ' num2str(f) 'of' 'Ori' num2str(o)]);

            for p=1:numPeriods

                spkb = convertSpikeTimes2Bins(spikeData(goodPos),analysisPeriodList{p},1000/Fs);
                lfpSpike = analogDataSpike(goodPos,xPos{p});
                lfp = analogData(goodPos,xPos{p});

                if removeEvokedResponse
                    lfpSpike =removeMeanResponse(lfpSpike);
                    lfp=removeMeanResponse(lfp);
                end

                % FFC and SFC
                [ffc(i,1,o,p,:),ffphi(i,1,o,p,:),~,~,~,ffcFreq]=coherencyc(lfpSpike',lfp',params); %#ok<*AGROW>
                [sfc(i,1,o,p,:),sfphi(i,1,o,p,:),~,~,~,sfcFreq]=coherencycpb(lfp',spkb,params);

                % Spike-triggered average
                [staVals(i,1,o,p,:),~,xsSTA] = getSTA(spikeData(goodPos),analogData(goodPos,:),analysisPeriodList{p},timeVals,[-0.1 0.1],1);

                % phase-amplitude coupling
                nFilterMethod=length(filterMethod);
%                 methods={'esc','mi','plv','tort','ozkurt'};
                nMethods=length(methodVar);
                for fi=1:nFilterMethod
                    filterName=filterMethod{fi};
                    if strcmp(filterName,'MP')
                  % Use the already filtered dat using Matching Pursuit
                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(spikeElectrode)],['elec' num2str(spikeElectrode) '.mat']));
                        analogDataSpikeMP=tmpMP.rSignal;
                        lfpSpikeMP = analogDataSpikeMP(goodPos,xPos{p});
                        

                        tmpMP = load(fullfile(folderSourceString,'ForToolbox', 'MP-master','MPdata', monkeyName,gridType,expDate,protocolName, ['elec' num2str(lfpElectrodes(i))],['elec' num2str(lfpElectrodes(i)) '.mat']));
                        analogDataMP=tmpMP.rSignal;
                        lfpMP = analogDataMP(goodPos,xPos{p});

                        if removeEvokedResponse
                            lfpSpikeMP =removeMeanResponse(lfpSpikeMP);
                            lfpMP=removeMeanResponse(lfpMP);
                        end

                        for m=1:nMethods
                            measure=methodVar{m};
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,1,o,p,fi,m,:,:), pval(i,1,o,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpikeMP', Fs, filterName, measure, ...
                                    lfpMP', ph_freq_vec, amp_freq_vec,  width, n_surr, epochframes, phFiltorder,ampFiltorder);
                        end
                    else
                        for m=1:nMethods
                            measure=methodVar{m};
                            % get PAC values and corresponding pvalues
                            % after permutation test
                            [pacmat(i,1,o,p,fi,m,:,:), pval(i,1,o,p,fi,m,:,:), freqvec_ph, freqvec_amp] = getPAC (lfpSpike', Fs, filterName, measure, ...
                                    lfp', ph_freq_vec, amp_freq_vec,  width, n_surr, epochframes, phFiltorder,ampFiltorder);
                        end
                    end
                end
                disp([measure ' for ' filterName 'for p' num2str(p)]);
            
            end
        end
    end
end
end

function y=removeMeanResponse(analogData)
y = analogData-repmat(mean(analogData),size(analogData,1),1);
end