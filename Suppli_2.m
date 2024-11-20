% Suppli 2: Plotting mean amplitude distribution of high frequency in each slow and fast gamma
% phase for D0 after MP case
% For both LFP and ECoG

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms
% Circular statistics: https://github.com/philippberens/circstat-matlab.git

clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Electrode choices
modality = 'LFP'; 

folderSave = 'savedData';

% Plot for Monkey 1 (M1)
% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol
% monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; 

% Plot Handles
phasePlotsM1 = getPlotHandles(2,1,[0.13 0.60 0.30 0.350],0.01,0.03);
plotHand = phasePlotsM1; 

% SF and Ori position for 'all-all' conditions
sfPos = 6; oriPos = 9;

% gamma range
gammaRange = {[40 72], [24 40]}; % Fast Gamma and Slow Gamma

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];
pacMethod = 'klmi'; 
filterName = 'fir'; 
nSurrogates = 0; % Number of surrogates
useMPFlag = 1; % Particular to MP filtered data
sVarName = 'sf';
electrodeDistanceVal = '0';

% distribution of amplitude providing frequency for each phase in phase
% signal
nBins=18; % as per Tort et al., 2010
position=zeros(1,nBins); 
winsize = 2*pi/nBins;
for nBin=1:nBins 
    position(nBin) = -pi+(nBin-1) * winsize; 
end

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% Frequency window used in calculating amplitude distribution
fGammaFreqWin = {[gammaRange{1}(1) gammaRange{1}(2) 157 487 ], [gammaRange{1}(1) gammaRange{1}(2) 87 157]};
sGammaFreqWin = {[gammaRange{2}(1) gammaRange{2}(2) 157 487 ], [gammaRange{2}(1) gammaRange{2}(2) 87 157]};

upperYlim =[0.05 0.065]; % High frequency >150Hz amplitude 
lowerYlim =[0.05 0.065]; % Low frequency <150Hz amplitude 
fGammaFreqWinColors = [153/255 102/255 204/255];  
sGammaFreqWinColors = [92/255 64/255 51/255];

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% load meanAmp
clear tmpData
tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' num2str(modality) '_d' num2str(electrodeDistanceVal) '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'meanAmp', 'centerAmpFreq', 'centerPhaseFreq');

%%%% Slow Gamma %%%%%
if mod(sGammaFreqWin{1,1}(1),4)==1
    xFreq1=sGammaFreqWin{1,1}(1)-1;
elseif mod(sGammaFreqWin{1,1}(1),4)==2
    xFreq1=sGammaFreqWin{1,1}(1)+2;
elseif mod(sGammaFreqWin{1,1}(1),4)==3
    xFreq1=sGammaFreqWin{1,1}(1)+1;
elseif mod(sGammaFreqWin{1,1}(1),4)==0
    xFreq1=sGammaFreqWin{1,1}(1);
end

if mod(sGammaFreqWin{1,1}(2),4)==1
    xFreq2=sGammaFreqWin{1,1}(2)-1;
elseif mod(sGammaFreqWin{1,1}(2),4)==2
    xFreq2=sGammaFreqWin{1,1}(2)+2;
elseif mod(sGammaFreqWin{1,1}(2),4)==3
    xFreq2=sGammaFreqWin{1,1}(2)+1;
elseif mod(sGammaFreqWin{1,1}(2),4)==0
    xFreq2=sGammaFreqWin{1,1}(2);
end
xFreqPosSG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

%%%%% Plotting Fast Gamma- Amp distribution %%%%%%
if mod(fGammaFreqWin{1,1}(1),4)==1
    xFreq1=fGammaFreqWin{1,1}(1)-1;
elseif mod(fGammaFreqWin{1,1}(1),4)==2
    xFreq1=fGammaFreqWin{1,1}(1)+2;
elseif mod(fGammaFreqWin{1,1}(1),4)==3
    xFreq1=fGammaFreqWin{1,1}(1)+1;
elseif mod(fGammaFreqWin{1,1}(1),4)==0
    xFreq1=fGammaFreqWin{1,1}(1);
end

if mod(fGammaFreqWin{1,1}(2),4)==1
    xFreq2=fGammaFreqWin{1,1}(2)-1;
elseif mod(fGammaFreqWin{1,1}(2),4)==2
    xFreq2=fGammaFreqWin{1,1}(2)+2;
elseif mod(fGammaFreqWin{1,1}(2),4)==3
    xFreq2=fGammaFreqWin{1,1}(2)+1;
elseif mod(fGammaFreqWin{1,1}(2),4)==0
    xFreq2=fGammaFreqWin{1,1}(2);
end
xFreqPosFG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

% > 150 Hz
yFreqPos1 = intersect(find(tmpData.centerAmpFreq>= sGammaFreqWin{1,1}(3)),find(tmpData.centerAmpFreq<= sGammaFreqWin{1,1}(4)));

clear normMeanAmpFG
ss=squeeze(tmpData.meanAmp(:,2,yFreqPos1,xFreqPosFG,:));
for nElec=1:size(ss,1)
    for iHF=1:size(ss,2)
        for iFG=1:size(ss,3)
             tmp2 = squeeze(ss(nElec,iHF,iFG,:));
             normMeanAmpFG(nElec,iHF,iFG,:) = tmp2/sum(tmp2);
        end
    end  
end

positionVal = position+pi;

% finding the angle at which the amplitude is peaked for each electrode
perElectrodeNormMeanAmpFG = squeeze(mean(normMeanAmpFG,[2 3]));
perElectrodePhaseForMaxValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
perElectrodeMaxValFGIdx = zeros(size(perElectrodeNormMeanAmpFG,1),1);
perElectrodeMaxValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
for nElec = 1:size(perElectrodeNormMeanAmpFG,1)
   perElectrodeMaxValFG(nElec,1) = max(perElectrodeNormMeanAmpFG(nElec,phasePos));
   perElectrodeMaxValFGIdx(nElec,1) = find(perElectrodeNormMeanAmpFG(nElec,phasePos) == perElectrodeMaxValFG(nElec,1));
   perElectrodePhaseForMaxValFG(nElec,1) = positionVal(perElectrodeMaxValFGIdx(nElec,1));
end

% significant test for each electrode's phase-mean ampli distribution
nBins = size(perElectrodeNormMeanAmpFG,2);

% initialise the array
klDivRawFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
pValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
tValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
for nElec = 1 : size(perElectrodeNormMeanAmpFG,1)
    klDivRawFG(nElec,1) = log(nBins) - (-sum(perElectrodeNormMeanAmpFG(nElec,:) .* log(perElectrodeNormMeanAmpFG(nElec,:))));

    % generate 100 surrogates for each phase-amplitude series
    nSurrogates = 1000;

    % randomise indices to generate shuffled data
    randomIdx = randi([1 nBins],nSurrogates,nBins);

    klDivSurr = zeros(nSurrogates,1);
    for nSurr = 1:nSurrogates
        clear shuffledData
        shuffledData = perElectrodeNormMeanAmpFG(nElec,randomIdx(nSurr,:));
        klDivSurr(nSurr,1) = log(nBins) - (-sum(shuffledData .* log(shuffledData)));
    end
    [~, pValue, ~, stats] = ttest2 (klDivSurr, klDivRawFG, 'Tail', 'both', 'Alpha', 0.01, 'Vartype', 'equal');
    tValFG(nElec,1) = stats.tstat;
    pValFG(nElec,1) = pValue;  
end

% plotting the phase-mean amplitude distribution for each electrode
cutOffFG = 1e-3;
for nElec=1:size(perElectrodeNormMeanAmpFG,1)
    if klDivRawFG(nElec,1) > cutOffFG
       plot(plotHand(1,1),position+pi,perElectrodeNormMeanAmpFG(nElec,phasePos),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
       hold(plotHand(1,1), "on");
       plot(plotHand(1,1), positionVal(perElectrodeMaxValFGIdx(nElec,1)), perElectrodeMaxValFG(nElec,1), 'Marker', 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)

    else
        plot(plotHand(1,1),position+pi,perElectrodeNormMeanAmpFG(nElec,phasePos),'Color',[0.85 0.85 0.85],'LineWidth',0.5)
    end
end
hold(plotHand(1,1), "on");

% Mean angle (circular mean)
avgPhaseForMaxValFG = (180/pi)*angle(mean(exp(1i*perElectrodePhaseForMaxValFG(find(klDivRawFG>cutOffFG),1))));

if (abs(avgPhaseForMaxValFG) ~= avgPhaseForMaxValFG)
    avgPhaseForMaxValFG = 360 + avgPhaseForMaxValFG;
end
disp(['Mean phase angle between FG and 150-500Hz in LFP M1=' num2str(avgPhaseForMaxValFG)]);
radAvgPhaseForMaxValFG = avgPhaseForMaxValFG*(pi/180);
radAvgPhaseForMaxValFGidx = find(positionVal<=radAvgPhaseForMaxValFG);
avgNormMeanAmpFG = squeeze(mean(normMeanAmpFG(find(klDivRawFG>cutOffFG),:,:,:),[1 2 3]));
avgNormMeanAmpFGtmp = avgNormMeanAmpFG(phasePos,1);
plot(plotHand(1,1),radAvgPhaseForMaxValFG, avgNormMeanAmpFGtmp(radAvgPhaseForMaxValFGidx(end)), 'Marker', 'o', 'Color', fGammaFreqWinColors, 'MarkerFaceColor', fGammaFreqWinColors, 'MarkerSize',8);
ylim(plotHand(1,1),upperYlim);
xlim(plotHand(1,1),[0 5.9341]);
set(plotHand(1,1),'xTick',[0 3.1416 5.9341 ]); set(plotHand(1,1),'xTickLabel',[]); 
box(plotHand(1,1),'off'); 
set(plotHand(1,1),'FontSize',12);
set(plotHand(1,1),'TickDir','out');
ylabel(plotHand(1,1),{' Amplitude of', '150 to 500 Hz'},'Fontsize',12);
text(4.0,0.064,['Fast Gamma (N=' num2str(size(perElectrodeNormMeanAmpFG(find(klDivRawFG>cutOffFG),:),1)) ')'],'Color', fGammaFreqWinColors, 'FontSize',12,'Parent',plotHand(1,1));

% < 150 Hz
yFreqPos2 = intersect(find(tmpData.centerAmpFreq>= sGammaFreqWin{1,2}(3)),find(tmpData.centerAmpFreq<= sGammaFreqWin{1,2}(4)));
clear normMeanAmpSG
ss=squeeze(tmpData.meanAmp(:,2,yFreqPos2,xFreqPosSG,:));
for nElec=1:size(ss,1)
    for iHF=1:size(ss,2)
        for iSG=1:size(ss,3)
             tmp1 = squeeze(ss(nElec,iHF,iSG,:));
             normMeanAmpSG(nElec,iHF,iSG,:) = tmp1/sum(tmp1);
        end
    end  
end

positionVal = position+pi;

% finding the angle at which the amplitude is peaked for each electrode
perElectrodeNormMeanAmpSG = squeeze(mean(normMeanAmpSG,[2 3]));
perElectrodePhaseForMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSGIdx = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
for nElec = 1:size(perElectrodeNormMeanAmpSG,1)
   perElectrodeMaxValSG(nElec,1) = max(perElectrodeNormMeanAmpSG(nElec,phasePos));
   perElectrodeMaxValSGIdx(nElec,1) = find(perElectrodeNormMeanAmpSG(nElec,phasePos) == perElectrodeMaxValSG(nElec,1));
   perElectrodePhaseForMaxValSG(nElec,1) = positionVal(perElectrodeMaxValSGIdx(nElec,1));
end

% significant test for each electrode's phase-mean ampli distribution
nBins = size(perElectrodeNormMeanAmpSG,2);

% initialise the array
klDivRawSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
pValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
tValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
for nElec = 1 : size(perElectrodeNormMeanAmpSG,1)
    klDivRawSG(nElec,1) = log(nBins) - (-sum(perElectrodeNormMeanAmpSG(nElec,:) .* log(perElectrodeNormMeanAmpSG(nElec,:))));

    % generate 100 surrogates for each phase-amplitude series
    nSurrogates = 1000;

    % randomise indices to generate shuffled data
    randomIdx = randi([1 nBins],nSurrogates,nBins);

    klDivSurr = zeros(nSurrogates,1);
    for nSurr = 1:nSurrogates
        clear shuffledData
        shuffledData = perElectrodeNormMeanAmpSG(nElec,randomIdx(nSurr,:));
        klDivSurr(nSurr,1) = log(nBins) - (-sum(shuffledData .* log(shuffledData)));
    end
     % Upaired t-test
      [~, pValue, ~, stats] = ttest2 (klDivSurr, klDivRawSG, 'Tail', 'both', 'Alpha', 0.01, 'Vartype', 'equal');
      tValSG(nElec,1) = stats.tstat;
      pValSG(nElec,1) = pValue;  
end

% plotting the phase-mean amplitude distribution for each electrode
cutOffSG = 1e-4;
for nElec=1:size(perElectrodeNormMeanAmpSG,1)
    if klDivRawSG(nElec,1) > cutOffSG
       plot(plotHand(2,1),position+pi,perElectrodeNormMeanAmpSG(nElec,phasePos),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
       hold(plotHand(2,1), "on");
       plot(plotHand(2,1), positionVal(perElectrodeMaxValSGIdx(nElec,1)), perElectrodeMaxValSG(nElec,1), 'Marker', 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)

    else
        plot(plotHand(2,1),position+pi,perElectrodeNormMeanAmpSG(nElec,phasePos),'Color',[0.85 0.85 0.85],'LineWidth',0.5)
    end
end
hold(plotHand(2,1), "on");

% Mean angle (circular mean)
avgPhaseForMaxValSG = (180/pi)*angle(mean(exp(1i*perElectrodePhaseForMaxValSG(find(klDivRawSG>cutOffSG),1))));

if (abs(avgPhaseForMaxValSG) ~= avgPhaseForMaxValSG)
    avgPhaseForMaxValSG = 360 + avgPhaseForMaxValSG;
end
disp(['Mean phase angle between SG and 80-150Hz in LFP M1=' num2str(avgPhaseForMaxValSG)]);
radAvgPhaseForMaxValSG = avgPhaseForMaxValSG*(pi/180);
radAvgPhaseForMaxValSGidx = find(positionVal<=radAvgPhaseForMaxValSG);
avgNormMeanAmpSG = squeeze(mean(normMeanAmpSG(find(klDivRawSG>cutOffSG),:,:,:),[1 2 3]));
avgNormMeanAmpSGtmp = avgNormMeanAmpSG(phasePos,1);
plot(plotHand(2,1),radAvgPhaseForMaxValSG, avgNormMeanAmpSGtmp(radAvgPhaseForMaxValSGidx(end)), 'Marker', 'o', 'Color', sGammaFreqWinColors, 'MarkerFaceColor', sGammaFreqWinColors, 'MarkerSize',8);
ylim(plotHand(2,1),upperYlim);
xlim(plotHand(2,1),[0 5.9341]);
set(plotHand(2,1),'xTick',[0 3.1416 5.9341 ]); set(plotHand(2,1),'xTickLabel',[]); 
set(plotHand(2,1),'xTickLabel',{'0','180','360'}); 
box(plotHand(2,1),'off'); 
set(plotHand(2,1),'FontSize',12);
set(plotHand(2,1),'TickDir','out');
ylabel(plotHand(2,1),{' Amplitude of', '150 to 500 Hz'},'Fontsize',12);
xlabel(plotHand(2,1),'Phase (degrees)','Fontsize',12);
text(4.0,0.064,['Slow Gamma(N=' num2str(size(perElectrodeNormMeanAmpSG(find(klDivRawSG>cutOffSG),:),1)) ')'],'Color', sGammaFreqWinColors, 'FontSize',12,'Parent',plotHand(2,1));

%************ Statistical Test *************
% Watson-William Test
[pvalSGFGM1, tableSGFGM1] = circ_wwtest(perElectrodePhaseForMaxValSG(find(klDivRawSG>cutOffSG),1), perElectrodePhaseForMaxValFG(find(klDivRawFG>cutOffFG),1));

annotation( 'textbox', 'String', 'Monkey 1', 'Color', 'black', ...
            'FontSize', 14,'FontWeight','Bold','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.26,0.3,0.9,0.7] )

%% Monkey 2
monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; 

% Plot Handles
phasePlotsM2 =getPlotHandles(2,1,[0.55 0.60 0.30 0.350],0.01,0.03);
plotHand = phasePlotsM2; 

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];
pacMethod = 'klmi'; 
filterName = 'fir'; 
nSurrogates = 0; % Number of surrogates
useMPFlag = 1; % Particular to MP filtered data
sVarName = 'sf';
electrodeDistanceVal = '0';

% distribution of amplitude providing frequency for each phase in phase
% signal
nBins=18; % as per Tort et al., 2010
position=zeros(1,nBins); 
winsize = 2*pi/nBins;
for nBin=1:nBins 
    position(nBin) = -pi+(nBin-1) * winsize; 
end

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% Frequency window used in calculating amplitude distribution
fGammaFreqWin = {[gammaRange{1}(1) gammaRange{1}(2) 157 487 ], [gammaRange{1}(1) gammaRange{1}(2) 87 157]};
sGammaFreqWin = {[gammaRange{2}(1) gammaRange{2}(2) 157 487 ], [gammaRange{2}(1) gammaRange{2}(2) 87 157]};

upperYlim =[0.05 0.065]; % High frequency >150Hz amplitude 
lowerYlim =[0.05 0.065]; % Low frequency <150Hz amplitude 
fGammaFreqWinColors = [153/255 102/255 204/255];  
sGammaFreqWinColors = [92/255 64/255 51/255];

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% load meanAmp
clear tmpData
tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' num2str(modality) '_d' num2str(electrodeDistanceVal) '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'meanAmp', 'centerAmpFreq', 'centerPhaseFreq');

%%%% Slow Gamma %%%%%
if mod(sGammaFreqWin{1,1}(1),4)==1
    xFreq1=sGammaFreqWin{1,1}(1)-1;
elseif mod(sGammaFreqWin{1,1}(1),4)==2
    xFreq1=sGammaFreqWin{1,1}(1)+2;
elseif mod(sGammaFreqWin{1,1}(1),4)==3
    xFreq1=sGammaFreqWin{1,1}(1)+1;
elseif mod(sGammaFreqWin{1,1}(1),4)==0
    xFreq1=sGammaFreqWin{1,1}(1);
end

if mod(sGammaFreqWin{1,1}(2),4)==1
    xFreq2=sGammaFreqWin{1,1}(2)-1;
elseif mod(sGammaFreqWin{1,1}(2),4)==2
    xFreq2=sGammaFreqWin{1,1}(2)+2;
elseif mod(sGammaFreqWin{1,1}(2),4)==3
    xFreq2=sGammaFreqWin{1,1}(2)+1;
elseif mod(sGammaFreqWin{1,1}(2),4)==0
    xFreq2=sGammaFreqWin{1,1}(2);
end
xFreqPosSG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

%%%%% Plotting Fast Gamma- Amp distribution %%%%%%
if mod(fGammaFreqWin{1,1}(1),4)==1
    xFreq1=fGammaFreqWin{1,1}(1)-1;
elseif mod(fGammaFreqWin{1,1}(1),4)==2
    xFreq1=fGammaFreqWin{1,1}(1)+2;
elseif mod(fGammaFreqWin{1,1}(1),4)==3
    xFreq1=fGammaFreqWin{1,1}(1)+1;
elseif mod(fGammaFreqWin{1,1}(1),4)==0
    xFreq1=fGammaFreqWin{1,1}(1);
end

if mod(fGammaFreqWin{1,1}(2),4)==1
    xFreq2=fGammaFreqWin{1,1}(2)-1;
elseif mod(fGammaFreqWin{1,1}(2),4)==2
    xFreq2=fGammaFreqWin{1,1}(2)+2;
elseif mod(fGammaFreqWin{1,1}(2),4)==3
    xFreq2=fGammaFreqWin{1,1}(2)+1;
elseif mod(fGammaFreqWin{1,1}(2),4)==0
    xFreq2=fGammaFreqWin{1,1}(2);
end
xFreqPosFG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

% > 150 Hz
yFreqPos1 = intersect(find(tmpData.centerAmpFreq>= sGammaFreqWin{1,1}(3)),find(tmpData.centerAmpFreq<= sGammaFreqWin{1,1}(4)));

clear normMeanAmpFG
ss=squeeze(tmpData.meanAmp(:,2,yFreqPos1,xFreqPosFG,:));
for nElec=1:size(ss,1)
    for iHF=1:size(ss,2)
        for iFG=1:size(ss,3)
             tmp2 = squeeze(ss(nElec,iHF,iFG,:));
             normMeanAmpFG(nElec,iHF,iFG,:) = tmp2/sum(tmp2);
        end
    end  
end

positionVal = position+pi;

% finding the angle at which the amplitude is peaked for each electrode
perElectrodeNormMeanAmpFG = squeeze(mean(normMeanAmpFG,[2 3]));
perElectrodePhaseForMaxValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
perElectrodeMaxValFGIdx = zeros(size(perElectrodeNormMeanAmpFG,1),1);
perElectrodeMaxValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
for nElec = 1:size(perElectrodeNormMeanAmpFG,1)
   perElectrodeMaxValFG(nElec,1) = max(perElectrodeNormMeanAmpFG(nElec,phasePos));
   perElectrodeMaxValFGIdx(nElec,1) = find(perElectrodeNormMeanAmpFG(nElec,phasePos) == perElectrodeMaxValFG(nElec,1));
   perElectrodePhaseForMaxValFG(nElec,1) = positionVal(perElectrodeMaxValFGIdx(nElec,1));
end

% significant test for each electrode's phase-mean ampli distribution
nBins = size(perElectrodeNormMeanAmpFG,2);
% initialise the array
klDivRawFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
pValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
tValFG = zeros(size(perElectrodeNormMeanAmpFG,1),1);
for nElec = 1 : size(perElectrodeNormMeanAmpFG,1)
    klDivRawFG(nElec,1) = log(nBins) - (-sum(perElectrodeNormMeanAmpFG(nElec,:) .* log(perElectrodeNormMeanAmpFG(nElec,:))));

    % generate 100 surrogates for each phase-amplitude series
    nSurrogates = 1000;

    % randomise indices to generate shuffled data
    randomIdx = randi([1 nBins],nSurrogates,nBins);

    klDivSurr = zeros(nSurrogates,1);
    for nSurr = 1:nSurrogates
        clear shuffledData
        shuffledData = perElectrodeNormMeanAmpFG(nElec,randomIdx(nSurr,:));
        klDivSurr(nSurr,1) = log(nBins) - (-sum(shuffledData .* log(shuffledData)));
    end

     % unpaired t-test 
      [~, pValue, ~, stats] = ttest2 (klDivSurr, klDivRawFG, 'Tail', 'both', 'Alpha', 0.01, 'Vartype', 'equal');
      tValFG(nElec,1) = stats.tstat;
      pValFG(nElec,1) = pValue;  
end

% plotting the phase-mean amplitude distribution for each electrode
cutOffFG = 0.75e-3;
for nElec=1:size(perElectrodeNormMeanAmpFG,1)
    if klDivRawFG(nElec,1) >= cutOffFG
       plot(plotHand(1,1),position+pi,perElectrodeNormMeanAmpFG(nElec,phasePos),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
       hold(plotHand(1,1), "on");
       plot(plotHand(1,1), positionVal(perElectrodeMaxValFGIdx(nElec,1)), perElectrodeMaxValFG(nElec,1), 'Marker', 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)

    elseif klDivRawFG(nElec,1) < cutOffFG
        plot(plotHand(1,1),position+pi,perElectrodeNormMeanAmpFG(nElec,phasePos),'Color',[0.85 0.85 0.85],'LineWidth',0.5)
    end
end
hold(plotHand(1,1), "on");

% Mean angle (circular mean)
avgPhaseForMaxValFG = (180/pi)*angle(mean(exp(1i*perElectrodePhaseForMaxValFG(find(klDivRawFG>cutOffFG),1))));

if (abs(avgPhaseForMaxValFG) ~= avgPhaseForMaxValFG)
    avgPhaseForMaxValFG = 360 + avgPhaseForMaxValFG;
end
disp(['Mean phase angle between FG and 150-500Hz in LFP M2=' num2str(avgPhaseForMaxValFG)]);
radAvgPhaseForMaxValFG = avgPhaseForMaxValFG*(pi/180);
radAvgPhaseForMaxValFGidx = find(positionVal<=radAvgPhaseForMaxValFG);
avgNormMeanAmpFG = squeeze(mean(normMeanAmpFG(find(klDivRawFG>cutOffFG),:,:,:),[1 2 3]));
avgNormMeanAmpFGtmp = avgNormMeanAmpFG(phasePos,1);
plot(plotHand(1,1),radAvgPhaseForMaxValFG, avgNormMeanAmpFGtmp(radAvgPhaseForMaxValFGidx(end)), 'Marker', 'o', 'Color', fGammaFreqWinColors, 'MarkerFaceColor', fGammaFreqWinColors, 'MarkerSize',8);
ylim(plotHand(1,1),upperYlim);
xlim(plotHand(1,1),[0 5.9341]);
set(plotHand(1,1),'xTick',[0 3.1416 5.9341 ]); set(plotHand(1,1),'xTickLabel',[]); 
box(plotHand(1,1),'off'); 
set(plotHand(1,1),'FontSize',12);
set(plotHand(1,1),'TickDir','out');
ylabel(plotHand(1,1),{' Amplitude of', '150 to 500 Hz'},'Fontsize',12);
text(4.0,0.064,['Fast Gamma (N=' num2str(size(perElectrodeNormMeanAmpFG(find(klDivRawFG>cutOffFG),:),1)) ')'],'Color', fGammaFreqWinColors, 'FontSize',12,'Parent',plotHand(1,1));

% < 150 Hz
yFreqPos2 = intersect(find(tmpData.centerAmpFreq>= sGammaFreqWin{1,2}(3)),find(tmpData.centerAmpFreq<= sGammaFreqWin{1,2}(4)));
clear normMeanAmpSG
ss=squeeze(tmpData.meanAmp(:,2,yFreqPos2,xFreqPosSG,:));
for nElec=1:size(ss,1)
    for iHF=1:size(ss,2)
        for iSG=1:size(ss,3)
             tmp1 = squeeze(ss(nElec,iHF,iSG,:));
             normMeanAmpSG(nElec,iHF,iSG,:) = tmp1/sum(tmp1);
        end
   end  
end

positionVal = position+pi;

% finding the angle at which the amplitude is peaked for each electrode
perElectrodeNormMeanAmpSG = squeeze(mean(normMeanAmpSG,[2 3]));
perElectrodePhaseForMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSGIdx = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
for nElec = 1:size(perElectrodeNormMeanAmpSG,1)
   perElectrodeMaxValSG(nElec,1) = max(perElectrodeNormMeanAmpSG(nElec,phasePos));
   perElectrodeMaxValSGIdx(nElec,1) = find(perElectrodeNormMeanAmpSG(nElec,phasePos) == perElectrodeMaxValSG(nElec,1));
   perElectrodePhaseForMaxValSG(nElec,1) = positionVal(perElectrodeMaxValSGIdx(nElec,1));
end

% significant test for each electrode's phase-mean ampli distribution
nBins = size(perElectrodeNormMeanAmpSG,2);

% initialise the array
klDivRawSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
pValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
tValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
for nElec = 1 : size(perElectrodeNormMeanAmpSG,1)
    klDivRawSG(nElec,1) = log(nBins) - (-sum(perElectrodeNormMeanAmpSG(nElec,:) .* log(perElectrodeNormMeanAmpSG(nElec,:))));

    % generate 100 surrogates for each phase-amplitude series
    nSurrogates = 1000;

    % randomise indices to generate shuffled data
    randomIdx = randi([1 nBins],nSurrogates,nBins);

    klDivSurr = zeros(nSurrogates,1);
    for nSurr = 1:nSurrogates
        clear shuffledData
        shuffledData = perElectrodeNormMeanAmpSG(nElec,randomIdx(nSurr,:));
        klDivSurr(nSurr,1) = log(nBins) - (-sum(shuffledData .* log(shuffledData)));
    end

     % unpaired t-test
      [~, pValue, ~, stats] = ttest2 (klDivSurr, klDivRawSG, 'Tail', 'both', 'Alpha', 0.01, 'Vartype', 'equal');
      tValSG(nElec,1) = stats.tstat;
      pValSG(nElec,1) = pValue;  
end

% plotting the phase-mean amplitude distribution for each electrode
cutOffSG = 1e-4;
for nElec=1:size(perElectrodeNormMeanAmpSG,1)
    if klDivRawSG(nElec,1) >= cutOffSG
       plot(plotHand(2,1),position+pi,perElectrodeNormMeanAmpSG(nElec,phasePos),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
       hold(plotHand(2,1), "on");
       plot(plotHand(2,1), positionVal(perElectrodeMaxValSGIdx(nElec,1)), perElectrodeMaxValSG(nElec,1), 'Marker', 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)

    elseif klDivRawSG(nElec,1) < cutOffSG
        plot(plotHand(2,1),position+pi,perElectrodeNormMeanAmpSG(nElec,phasePos),'Color',[0.85 0.85 0.85],'LineWidth',0.5)
    end
end
hold(plotHand(2,1), "on");

% Mean angle (circular mean)
avgPhaseForMaxValSG = (180/pi)*angle(mean(exp(1i*perElectrodePhaseForMaxValSG(find(klDivRawSG>cutOffSG),1))));

if (abs(avgPhaseForMaxValSG) ~= avgPhaseForMaxValSG)
    avgPhaseForMaxValSG = 360 + avgPhaseForMaxValSG;
end
disp(['Mean phase angle between SG and 80-150Hz in LFP M2=' num2str(avgPhaseForMaxValSG)]);
radAvgPhaseForMaxValSG = avgPhaseForMaxValSG*(pi/180);
radAvgPhaseForMaxValSGidx = find(positionVal<=radAvgPhaseForMaxValSG);
avgNormMeanAmpSG = squeeze(mean(normMeanAmpSG(find(klDivRawSG>cutOffSG),:,:,:),[1 2 3]));
avgNormMeanAmpSGtmp = avgNormMeanAmpSG(phasePos,1);
plot(plotHand(2,1),radAvgPhaseForMaxValSG, avgNormMeanAmpSGtmp(radAvgPhaseForMaxValSGidx(end)), 'Marker', 'o', 'Color', sGammaFreqWinColors, 'MarkerFaceColor', sGammaFreqWinColors, 'MarkerSize',8);
% plot(plotHand(2,1),[1 1]*positionVal(radAvgPhaseForMaxValSGidx(end)), [0 1]*avgNormMeanAmpSGtmp(radAvgPhaseForMaxValSGidx(end)), 'LineStyle', '--', 'Color', sGammaFreqWinColors);                          % Vertical Dashed Line
ylim(plotHand(2,1),upperYlim);
xlim(plotHand(2,1),[0 5.9341]);
set(plotHand(2,1),'xTick',[0 3.1416 5.9341 ]); set(plotHand(2,1),'xTickLabel',[]); 
set(plotHand(2,1),'xTickLabel',{'0','180','360'}); 
box(plotHand(2,1),'off'); 
set(plotHand(2,1),'FontSize',12);
set(plotHand(2,1),'TickDir','out');
ylabel(plotHand(2,1),{' Amplitude of', '150 to 500 Hz'},'Fontsize',12);
xlabel(plotHand(2,1),'Phase (degrees)','Fontsize',12);
text(4.0,0.064,['Slow Gamma(N=' num2str(size(perElectrodeNormMeanAmpSG(find(klDivRawSG>cutOffSG),:),1)) ')'],'Color', sGammaFreqWinColors, 'FontSize',12,'Parent',plotHand(2,1));

%************ Statistical Test *************
% Watson-William Test
[pvalSGFGM2, tableSGFGM2] = circ_wwtest(perElectrodePhaseForMaxValSG(find(klDivRawSG>cutOffSG),1), perElectrodePhaseForMaxValFG(find(klDivRawFG>cutOffFG),1));

annotation( 'textbox', 'String', 'Monkey 2', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight','Bold','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.68,0.3,0.9,0.7] )
annotation( 'textbox', 'String', 'LFP', 'Color', 'black', ...
            'FontSize', 16, 'FontAngle','italic','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.44,0.3,0.9,0.5] )

%%
%%%%%%%%%%%%%% ECoG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrode choices
modality = 'ECoG'; %'LFP' or 'ECoG'

%% Monkey 1
folderSave = 'savedData';

% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol

phasePlotsM1 = getPlotHandles(1,1,[0.13 0.32 0.30 0.18],0.01,0.03);
plotHand = phasePlotsM1; 

% SF and Ori position for 'all-all' conditions
sfPos = 6; oriPos = 9;

% gamma range
gammaRange = {[40 72], [24 40]}; % Fast Gamma and Slow Gamma

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];
pacMethod = 'klmi'; 
filterName = 'fir'; 
nSurrogates = 0; % Number of surrogates
useMPFlag = 1; % Particular to MP filtered data
sVarName = 'sf';
electrodeDistanceVal = '0';

% distribution of amplitude providing frequency for each phase in phase
% signal
nBins=18; % as per Tort et al., 2010
position=zeros(1,nBins); 
winsize = 2*pi/nBins;
for nBin=1:nBins 
    position(nBin) = -pi+(nBin-1) * winsize; 
end

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% Frequency window used in calculating amplitude distribution
fGammaFreqWin = {[gammaRange{1}(1) gammaRange{1}(2) 157 487 ], [gammaRange{1}(1) gammaRange{1}(2) 87 157]};
sGammaFreqWin = {[gammaRange{2}(1) gammaRange{2}(2) 157 487 ], [gammaRange{2}(1) gammaRange{2}(2) 87 157]};

upperYlim =[0.045 0.07]; % High frequency >150Hz amplitude 
lowerYlim =[0.045 0.07]; % Low frequency <150Hz amplitude 
fGammaFreqWinColors = [153/255 102/255 204/255];  
sGammaFreqWinColors = [92/255 64/255 51/255];

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% load meanAmp
clear tmpData
tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' num2str(modality) '_d' num2str(electrodeDistanceVal) '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'meanAmp', 'centerAmpFreq', 'centerPhaseFreq');

%%%% Slow Gamma %%%%%
if mod(sGammaFreqWin{1,1}(1),4)==1
    xFreq1=sGammaFreqWin{1,1}(1)-1;
elseif mod(sGammaFreqWin{1,1}(1),4)==2
    xFreq1=sGammaFreqWin{1,1}(1)+2;
elseif mod(sGammaFreqWin{1,1}(1),4)==3
    xFreq1=sGammaFreqWin{1,1}(1)+1;
elseif mod(sGammaFreqWin{1,1}(1),4)==0
    xFreq1=sGammaFreqWin{1,1}(1);
end

if mod(sGammaFreqWin{1,1}(2),4)==1
    xFreq2=sGammaFreqWin{1,1}(2)-1;
elseif mod(sGammaFreqWin{1,1}(2),4)==2
    xFreq2=sGammaFreqWin{1,1}(2)+2;
elseif mod(sGammaFreqWin{1,1}(2),4)==3
    xFreq2=sGammaFreqWin{1,1}(2)+1;
elseif mod(sGammaFreqWin{1,1}(2),4)==0
    xFreq2=sGammaFreqWin{1,1}(2);
end
xFreqPosSG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

%%%%% Plotting Fast Gamma- Amp distribution %%%%%%
if mod(fGammaFreqWin{1,1}(1),4)==1
    xFreq1=fGammaFreqWin{1,1}(1)-1;
elseif mod(fGammaFreqWin{1,1}(1),4)==2
    xFreq1=fGammaFreqWin{1,1}(1)+2;
elseif mod(fGammaFreqWin{1,1}(1),4)==3
    xFreq1=fGammaFreqWin{1,1}(1)+1;
elseif mod(fGammaFreqWin{1,1}(1),4)==0
    xFreq1=fGammaFreqWin{1,1}(1);
end

if mod(fGammaFreqWin{1,1}(2),4)==1
    xFreq2=fGammaFreqWin{1,1}(2)-1;
elseif mod(fGammaFreqWin{1,1}(2),4)==2
    xFreq2=fGammaFreqWin{1,1}(2)+2;
elseif mod(fGammaFreqWin{1,1}(2),4)==3
    xFreq2=fGammaFreqWin{1,1}(2)+1;
elseif mod(fGammaFreqWin{1,1}(2),4)==0
    xFreq2=fGammaFreqWin{1,1}(2);
end
xFreqPosFG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

 % < 150 Hz
yFreqPos2 = intersect(find(tmpData.centerAmpFreq>= sGammaFreqWin{1,2}(3)),find(tmpData.centerAmpFreq<= sGammaFreqWin{1,2}(4)));
clear normMeanAmpSG
ss=squeeze(tmpData.meanAmp(:,2,yFreqPos2,xFreqPosSG,:));
for nElec=1:size(ss,1)
    for iHF=1:size(ss,2)
        for iSG=1:size(ss,3)
             tmp1 = squeeze(ss(nElec,iHF,iSG,:));
             normMeanAmpSG(nElec,iHF,iSG,:) = tmp1/sum(tmp1);
        end
    end  
end

positionVal = position+pi;

% finding the angle at which the amplitude is peaked for each electrode
perElectrodeNormMeanAmpSG = squeeze(mean(normMeanAmpSG,[2 3]));
perElectrodePhaseForMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSGIdx = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
for nElec = 1:size(perElectrodeNormMeanAmpSG,1)
   perElectrodeMaxValSG(nElec,1) = max(perElectrodeNormMeanAmpSG(nElec,phasePos));
   perElectrodeMaxValSGIdx(nElec,1) = find(perElectrodeNormMeanAmpSG(nElec,phasePos) == perElectrodeMaxValSG(nElec,1));
   perElectrodePhaseForMaxValSG(nElec,1) = positionVal(perElectrodeMaxValSGIdx(nElec,1));
end

% Mean angle (circular mean)
avgPhaseForMaxValSG = (180/pi)*angle(mean(exp(1i*perElectrodePhaseForMaxValSG)));

% plotting the phase-mean amplitude distribution for each electrode
for nElec=1:size(perElectrodeNormMeanAmpSG,1)
    plot(plotHand(1,1),position+pi,perElectrodeNormMeanAmpSG(nElec,phasePos),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
    hold(plotHand(1,1), "on");
    plot(plotHand(1,1), positionVal(perElectrodeMaxValSGIdx(nElec,1)), perElectrodeMaxValSG(nElec,1), 'Marker', 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
end
hold(plotHand(1,1), "on");
avgNormMeanAmpSG = squeeze(mean(normMeanAmpSG,[1 2 3]));
if (abs(avgPhaseForMaxValSG) ~= avgPhaseForMaxValSG)
    avgPhaseForMaxValSG = 360 + avgPhaseForMaxValSG;
end
disp(['Mean phase angle between SG and 80-150 Hz in ECoG M1=' num2str(avgPhaseForMaxValSG)]);

radAvgPhaseForMaxValSG = avgPhaseForMaxValSG*(pi/180);
radAvgPhaseForMaxValSGidx = find(positionVal<=radAvgPhaseForMaxValSG);
avgNormMeanAmpSGtmp = avgNormMeanAmpSG(phasePos,1);
plot(plotHand(1,1),positionVal(radAvgPhaseForMaxValSGidx(end)), avgNormMeanAmpSGtmp(radAvgPhaseForMaxValSGidx(end)), 'Marker', 'o', 'Color', sGammaFreqWinColors, 'MarkerFaceColor', sGammaFreqWinColors, 'MarkerSize',8);

ylim(plotHand(1,1),upperYlim);
xlim(plotHand(1,1),[0 5.9341]);
set(plotHand(1,1),'xTick',[0 3.1416 5.9341 ]); set(plotHand(1,1),'xTickLabel',[]); 
set(plotHand(1,1),'xTickLabel',{'0','180','360'});
box(plotHand(1,1),'off'); 
set(plotHand(1,1),'FontSize',12);
set(plotHand(1,1),'TickDir','out');
ylabel(plotHand(1,1),{' Amplitude of', '80 to 150 Hz'},'Fontsize',12);
xlabel(plotHand(1,1),'Phase (degrees)','Fontsize',12);
text(4.0,0.069,['Slow Gamma(N=' num2str(size(perElectrodeNormMeanAmpSG,1)) ')'],'Color', sGammaFreqWinColors, 'FontSize',12,'Parent',plotHand(1,1));

%% Monkey 2

% Protocol details
monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; 

% Plot Handles
phasePlotsM2 =getPlotHandles(1,1,[0.55 0.32 0.30 0.18],0.01,0.03);
plotHand = phasePlotsM2; 

% distribution of amplitude providing frequency for each phase in phase
% signal
nBins=18; % as per Tort et al., 2010
position=zeros(1,nBins); 
winsize = 2*pi/nBins;
for nBin=1:nBins 
    position(nBin) = -pi+(nBin-1) * winsize; 
end

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% Frequency window used in calculating amplitude distribution
fGammaFreqWin = {[gammaRange{1}(1) gammaRange{1}(2) 157 487 ], [gammaRange{1}(1) gammaRange{1}(2) 87 157]};
sGammaFreqWin = {[gammaRange{2}(1) gammaRange{2}(2) 157 487 ], [gammaRange{2}(1) gammaRange{2}(2) 87 157]};

upperYlim =[0.045 0.065]; % High frequency >150Hz amplitude 
lowerYlim =[0.045 0.065]; % Low frequency <150Hz amplitude 
fGammaFreqWinColors = [153/255 102/255 204/255];  
sGammaFreqWinColors = [92/255 64/255 51/255];

phasePos1 = intersect(find(position>= 0),find(position<=  pi));
phasePos2 = intersect(find(position< 0),find(position>=  -pi));
phasePos = [phasePos1 phasePos2];

% load meanAmp
clear tmpData
tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' num2str(modality) '_d' num2str(electrodeDistanceVal) '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'meanAmp', 'centerAmpFreq', 'centerPhaseFreq');
%%%% Slow Gamma %%%%%
if mod(sGammaFreqWin{1,1}(1),4)==1
    xFreq1=sGammaFreqWin{1,1}(1)-1;
elseif mod(sGammaFreqWin{1,1}(1),4)==2
    xFreq1=sGammaFreqWin{1,1}(1)+2;
elseif mod(sGammaFreqWin{1,1}(1),4)==3
    xFreq1=sGammaFreqWin{1,1}(1)+1;
elseif mod(sGammaFreqWin{1,1}(1),4)==0
    xFreq1=sGammaFreqWin{1,1}(1);
end

if mod(sGammaFreqWin{1,1}(2),4)==1
    xFreq2=sGammaFreqWin{1,1}(2)-1;
elseif mod(sGammaFreqWin{1,1}(2),4)==2
    xFreq2=sGammaFreqWin{1,1}(2)+2;
elseif mod(sGammaFreqWin{1,1}(2),4)==3
    xFreq2=sGammaFreqWin{1,1}(2)+1;
elseif mod(sGammaFreqWin{1,1}(2),4)==0
    xFreq2=sGammaFreqWin{1,1}(2);
end
xFreqPosSG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

%%%%% Plotting Fast Gamma- Amp distribution %%%%%%
if mod(fGammaFreqWin{1,1}(1),4)==1
    xFreq1=fGammaFreqWin{1,1}(1)-1;
elseif mod(fGammaFreqWin{1,1}(1),4)==2
    xFreq1=fGammaFreqWin{1,1}(1)+2;
elseif mod(fGammaFreqWin{1,1}(1),4)==3
    xFreq1=fGammaFreqWin{1,1}(1)+1;
elseif mod(fGammaFreqWin{1,1}(1),4)==0
    xFreq1=fGammaFreqWin{1,1}(1);
end

if mod(fGammaFreqWin{1,1}(2),4)==1
    xFreq2=fGammaFreqWin{1,1}(2)-1;
elseif mod(fGammaFreqWin{1,1}(2),4)==2
    xFreq2=fGammaFreqWin{1,1}(2)+2;
elseif mod(fGammaFreqWin{1,1}(2),4)==3
    xFreq2=fGammaFreqWin{1,1}(2)+1;
elseif mod(fGammaFreqWin{1,1}(2),4)==0
    xFreq2=fGammaFreqWin{1,1}(2);
end
xFreqPosFG = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));

 % < 150 Hz
yFreqPos2 = intersect(find(tmpData.centerAmpFreq>= sGammaFreqWin{1,2}(3)),find(tmpData.centerAmpFreq<= sGammaFreqWin{1,2}(4)));
clear normMeanAmpSG
ss=squeeze(tmpData.meanAmp(:,2,yFreqPos2,xFreqPosSG,:));
for nElec=1:size(ss,1)
    for iHF=1:size(ss,2)
        for iSG=1:size(ss,3)
             tmp1 = squeeze(ss(nElec,iHF,iSG,:));
             normMeanAmpSG(nElec,iHF,iSG,:) = tmp1/sum(tmp1);
        end
    end  
end

positionVal = position+pi;

% finding the angle at which the amplitude is peaked for each electrode
perElectrodeNormMeanAmpSG = squeeze(mean(normMeanAmpSG,[2 3]));
perElectrodePhaseForMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSGIdx = zeros(size(perElectrodeNormMeanAmpSG,1),1);
perElectrodeMaxValSG = zeros(size(perElectrodeNormMeanAmpSG,1),1);
for nElec = 1:size(perElectrodeNormMeanAmpSG,1)
   perElectrodeMaxValSG(nElec,1) = max(perElectrodeNormMeanAmpSG(nElec,phasePos));
   perElectrodeMaxValSGIdx(nElec,1) = find(perElectrodeNormMeanAmpSG(nElec,phasePos) == perElectrodeMaxValSG(nElec,1));
   perElectrodePhaseForMaxValSG(nElec,1) = positionVal(perElectrodeMaxValSGIdx(nElec,1));
end

% Mean angle (circular mean)
avgPhaseForMaxValSG = (180/pi)*angle(mean(exp(1i*perElectrodePhaseForMaxValSG)));

% plotting the phase-mean amplitude distribution for each electrode
for nElec=1:size(perElectrodeNormMeanAmpSG,1)
    plot(plotHand(1,1),position+pi,perElectrodeNormMeanAmpSG(nElec,phasePos),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
    hold(plotHand(1,1), "on");
    plot(plotHand(1,1), positionVal(perElectrodeMaxValSGIdx(nElec,1)), perElectrodeMaxValSG(nElec,1), 'Marker', 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
end
hold(plotHand(1,1), "on");
avgNormMeanAmpSG = squeeze(mean(normMeanAmpSG,[1 2 3]));
if (abs(avgPhaseForMaxValSG) ~= avgPhaseForMaxValSG)
    avgPhaseForMaxValSG = 360 + avgPhaseForMaxValSG;
end
disp(['Mean phase angle between SG and 80-150 Hz in ECoG M2=' num2str(avgPhaseForMaxValSG)]);

radAvgPhaseForMaxValSG = avgPhaseForMaxValSG*(pi/180);
radAvgPhaseForMaxValSGidx = find(positionVal<=radAvgPhaseForMaxValSG);
avgNormMeanAmpSGtmp = avgNormMeanAmpSG(phasePos,1);
plot(plotHand(1,1),positionVal(radAvgPhaseForMaxValSGidx(end)), avgNormMeanAmpSGtmp(radAvgPhaseForMaxValSGidx(end)), 'Marker', 'o', 'Color', sGammaFreqWinColors, 'MarkerFaceColor', sGammaFreqWinColors, 'MarkerSize',8);

ylim(plotHand(1,1),upperYlim);
xlim(plotHand(1,1),[0 5.9341]);
set(plotHand(1,1),'xTick',[0 3.1416 5.9341 ]); set(plotHand(1,1),'xTickLabel',[]); 
set(plotHand(1,1),'xTickLabel',{'0','180','360'});
box(plotHand(1,1),'off'); 
set(plotHand(1,1),'FontSize',12);
set(plotHand(1,1),'TickDir','out');
ylabel(plotHand(1,1),{' Amplitude of', '80 to 150 Hz'},'Fontsize',12);
xlabel(plotHand(1,1),'Phase (degrees)','Fontsize',12);
text(4.0,0.064,['Slow Gamma(N=' num2str(size(perElectrodeNormMeanAmpSG,1)) ')'],'Color', sGammaFreqWinColors, 'FontSize',12,'Parent',plotHand(1,1));

annotation( 'textbox', 'String', 'ECoG', 'Color', 'black', ...
            'FontSize', 16, 'FontAngle','italic','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.44,0.33,0.9,0.1] )
