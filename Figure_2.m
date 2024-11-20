% Figure 2: Plotting for Spike Triggered average (STA) and Spike Field 
% Coherence (SFC) for all-all conditions for Local Field Potential (LFP), under three cases
% 1. D0 (single electrode LFP)
% 2. D0 after MP (single electrode LFP after Matching Pursuit(MP))
% 3. D400 (inter-electrode within 400 micrometer distance)

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrode choices
modality = 'LFP'; 

folderSave = 'savedData';

% Color vector for three curves (D0, D0 after MP, D400)
colorVals = [1 0 1; 0.2 0.8 0.5; 1 0.6824 0.3176];

% SF and Ori position for 'all-all' conditions
sfPos = 6; oriPos = 9;

% Fast Gamma and Slow Gamma range, respectively
gammaRange = {[40 70], [25 40]}; % Chosen as per the peaks visible in PSD in both the Monkeys

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];
pacMethod = 'klmi'; 
filterName = 'fir'; 
nSurrogates = 0; % Number of surrogates
sVarName = 'sf';

%% Monkey 1

% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol

% Plot Handles
psdPlotsM1 = getPlotHandles(2,1,[0.09 0.38 0.35 0.55],0.04,0.09);
plotHand = psdPlotsM1;

% Load STA and SFC values for D0 
useMPFlag = '0'; 
electrodeDistanceVal = '0';
tmpData1 = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'staVals', 'xsSTA', 'sfc', 'sfcFreq');

% Load STA and SFC values for D0 after MP
useMPFlag = '1'; 
electrodeDistanceVal = '0';
tmpData2 = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'staVals', 'xsSTA', 'sfc', 'sfcFreq');

% Load STA and SFC for D400 
useMPFlag = '0'; 
electrodeDistanceVal = '400';
tmpData3 = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'staVals', 'xsSTA', 'sfc', 'sfcFreq');

%%%%%%%%%%% Plotting STA for D0, D0 after MP and D400 %%%%%%%%%%%%%%%%%%%%%
plot(plotHand(1,1),1000*tmpData1.xsSTA,squeeze(mean(tmpData1.staVals(:,2,:),1)),'color',colorVals(1,:),'LineWidth',1.5);
hold(plotHand(1,1), "on");
plot(plotHand(1,1),1000*tmpData2.xsSTA,squeeze(mean(tmpData2.staVals(:,2,:),1)),'color',colorVals(2,:),'LineWidth',1.5);
plot(plotHand(1,1),1000*tmpData3.xsSTA,squeeze(mean(tmpData3.staVals(:,2,:),1)),'color',colorVals(3,:),'LineWidth',1.5);
set(plotHand(1,1),'NextPlot','add'); axis(plotHand(1,1),[tmpData1.xsSTA(1)*10^3 tmpData1.xsSTA(end)*10^3 -40 0]); set(plotHand(1,1),'TickDir','out'); box(plotHand(1,1),'off');

%%%%%%%%%%% Plotting SFC for D0, D0 after MP and D400 %%%%%%%%%%%%%%%%%%%%%
plot(plotHand(2,1),tmpData1.sfcFreq,squeeze(mean(tmpData1.sfc(:,2,:),1)),'color',colorVals(1,:),'LineWidth',1.5);
hold(plotHand(2,1), "on");
plot(plotHand(2,1),tmpData2.sfcFreq,squeeze(mean(tmpData2.sfc(:,2,:),1)),'color',colorVals(2,:),'LineWidth',1.5);
plot(plotHand(2,1),tmpData3.sfcFreq,squeeze(mean(tmpData3.sfc(:,2,:),1)),'color',colorVals(3,:),'LineWidth',1.5);
set(plotHand(2,1),'NextPlot','add'); axis(plotHand(2,1),[0 250 0 0.2]); set(plotHand(2,1),'TickDir','out'); box(plotHand(2,1),'off');
set(plotHand(2,1),'Xscale','log');
set(plotHand(1,1),'xTick',[-80 -40 0 40 80]); set(plotHand(1,1),'xTickLabel',[-80 -40 0 40 80]);
set(plotHand(2,1),'xTick',[10 30 70 150]); set(plotHand(2,1),'xTickLabel',[10 30 70 150]);
set(plotHand(1,1),'yTick',[-40 -20 0]); set(plotHand(1,1),'yTickLabel',[-40 -20 0]);
set(plotHand(2,1),'yTick',[0 0.1 0.2]); set(plotHand(2,1),'yTickLabel',[0 0.1 0.2]);
set(plotHand(1,1),'FontSize',15);
set(plotHand(2,1),'FontSize',15);
set(plotHand(1,1),'LineWidth',1);
set(plotHand(2,1),'LineWidth',1);
% Fast gamma
xline(plotHand(2,1),gammaRange{1}(1),'Color',[153/255 102/255 204/255],'LineWidth',1);xline(plotHand(2,1),gammaRange{1}(2),'Color',[153/255 102/255 204/255],'LineWidth',1);
% slow gamma
xline(plotHand(2,1),gammaRange{2}(1),'Color',[92/255 64/255 51/255],'LineWidth',1);xline(plotHand(2,1),gammaRange{2}(2),'Color',[92/255 64/255 0/51],'LineWidth',1);
ylabel(plotHand(1,1), 'STA (\muV)','FontSize',15); xlabel(plotHand(1,1), 'Time (ms)','Fontsize',15); 
ylabel(plotHand(2,1), 'SFC','FontSize',15); xlabel(plotHand(2,1), 'Frequency (Hz)','Fontsize',15);
text(40,-23,['D_{0}(N=' num2str(size(tmpData1.sfc,1)) ')'],'Color', colorVals(1,:), 'FontSize',13,'Parent',plotHand(1,1));
text(40,-28,['D_{0} after MP(N=' num2str(size(tmpData2.sfc,1)) ')'],'Color', colorVals(2,:), 'FontSize',13,'Parent',plotHand(1,1));
text(40,-33,['D_{400}(N=' num2str(size(tmpData3.sfc,1)) ')'],'Color', colorVals(3,:), 'FontSize',13,'Parent',plotHand(1,1));

annotation( 'textbox', 'String', 'Monkey 1', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight', 'Bold', 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.24,0.9,0.9,0.09], 'rotation',0 )

%% Monkey 2

% Protocol details
monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; 

% Plot Handles
psdPlotsM2 = getPlotHandles(2,1,[0.50 0.38 0.35 0.55],0.04,0.09);
plotHand = psdPlotsM2;

% Load STA  and SFC values for D0 
useMPFlag = '0'; 
electrodeDistanceVal = '0';
tmpData1 = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'staVals', 'xsSTA', 'sfc', 'sfcFreq');

% Load STA  and SFC values for D0 after MP
useMPFlag = '1'; 
electrodeDistanceVal = '0';
tmpData2 = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'staVals', 'xsSTA', 'sfc', 'sfcFreq');

% Load STA and SFC for D400 
useMPFlag = '0'; 
electrodeDistanceVal = '400';
tmpData3 = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'staVals', 'xsSTA', 'sfc', 'sfcFreq');

%%%%%%%%%%% Plotting STA for D0, D0 after MP and D400 %%%%%%%%%%%%%%%%%%%%%
plot(plotHand(1,1),1000*tmpData1.xsSTA,squeeze(mean(tmpData1.staVals(:,2,:),1)),'color',colorVals(1,:),'LineWidth',1.5);
hold(plotHand(1,1), "on");
plot(plotHand(1,1),1000*tmpData2.xsSTA,squeeze(mean(tmpData2.staVals(:,2,:),1)),'color',colorVals(2,:),'LineWidth',1.5);
plot(plotHand(1,1),1000*tmpData3.xsSTA,squeeze(mean(tmpData3.staVals(:,2,:),1)),'color',colorVals(3,:),'LineWidth',1.5);
set(plotHand(1,1),'NextPlot','add'); axis(plotHand(1,1),[tmpData1.xsSTA(1)*10^3 tmpData1.xsSTA(end)*10^3 -20 10]); set(plotHand(1,1),'TickDir','out'); box(plotHand(1,1),'off');

%%%%%%%%%%% Plotting SFC for D0, D0 after MP and D400 %%%%%%%%%%%%%%%%%%%%%
plot(plotHand(2,1),tmpData1.sfcFreq,squeeze(mean(tmpData1.sfc(:,2,:),1)),'color',colorVals(1,:),'LineWidth',1.5);
hold(plotHand(2,1), "on");
plot(plotHand(2,1),tmpData2.sfcFreq,squeeze(mean(tmpData2.sfc(:,2,:),1)),'color',colorVals(2,:),'LineWidth',1.5);
plot(plotHand(2,1),tmpData3.sfcFreq,squeeze(mean(tmpData3.sfc(:,2,:),1)),'color',colorVals(3,:),'LineWidth',1.5);
set(plotHand(2,1),'NextPlot','add'); axis(plotHand(2,1),[0 250 0 0.2]); set(plotHand(2,1),'TickDir','out'); box(plotHand(2,1),'off');
set(plotHand(2,1),'Xscale','log');
set(plotHand(1,1),'xTick',[-80 -40 0 40 80]); set(plotHand(1,1),'xTickLabel',[-80 -40 0 40 80]);
set(plotHand(2,1),'xTick',[10 30 70 150]); set(plotHand(2,1),'xTickLabel',[10 30 70 150]);
set(plotHand(2,1),'yTick',[0 0.1 0.2]); set(plotHand(2,1),'yTickLabel',[0 0.1 0.2]);
set(plotHand(1,1),'FontSize',15);
set(plotHand(2,1),'FontSize',15);
set(plotHand(1,1),'LineWidth',1);
set(plotHand(2,1),'LineWidth',1);
% Fast gamma
xline(plotHand(2,1),gammaRange{1}(1),'Color',[153/255 102/255 204/255],'LineWidth',1);xline(plotHand(2,1),gammaRange{1}(2),'Color',[153/255 102/255 204/255],'LineWidth',1);
% Slow gamma
xline(plotHand(2,1),gammaRange{2}(1),'Color',[92/255 64/255 51/255],'LineWidth',1);xline(plotHand(2,1),gammaRange{2}(2),'Color',[92/255 64/255 51/255],'LineWidth',1);
ylabel(plotHand(1,1), 'STA (\muV)','FontSize',15); xlabel(plotHand(1,1), 'Time (ms)','Fontsize',15); 
ylabel(plotHand(2,1), 'SFC','FontSize',15); xlabel(plotHand(2,1), 'Frequency (Hz)','Fontsize',15);
text(40,-8,['N=' num2str(size(tmpData1.sfc,1))],'Color', colorVals(1,:), 'FontSize',13,'Parent',plotHand(1,1));
text(40,-12,['N=' num2str(size(tmpData2.sfc,1))],'Color', colorVals(2,:), 'FontSize',13,'Parent',plotHand(1,1));
text(40,-16,['N=' num2str(size(tmpData3.sfc,1))],'Color', colorVals(3,:), 'FontSize',13,'Parent',plotHand(1,1));

annotation( 'textbox', 'String', 'Monkey 2', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight', 'Bold', 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.65,0.9,0.9,0.09], 'rotation',0 )