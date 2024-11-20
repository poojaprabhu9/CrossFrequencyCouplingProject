% Figure 1: Plotting for all-all conditions
% Plotting the PSD values for Baseline and Stimulus period and delta PSD
% for LFP and ECoG

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms

clear;
close all

%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrode choices
modality = 'LFP'; %'LFP' or 'ECoG'

% Signal processing choices
tapers = [2 3];

folderSave = 'savedData';

% SF and Ori position for 'all-all' conditions
sfPos = 6; oriPos = 9;

% Fast Gamma and Slow Gamma range, respectively
gammaRange = {[40 70], [25 40]}; % Fast Gamma and Slow Gamma; Chosen as per the peaks visible in PSD in both the Monkeys

%% Monkey 1

% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol

% Plot Handles
psdPlotsM1 = getPlotHandles(2,1,[0.09 0.34 0.15 0.6],0.01,0.03);
deltaPsdPlotsM1 = getPlotHandles(1,1,[0.3 0.34 0.15 0.6],0.01,0.01);
plotHand = psdPlotsM1;
deltaPlotHand = deltaPsdPlotsM1;

% load PSD data
psdData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_spike' modality '_PSD_For_All_Conditions.mat']));

%%%%%%%%%%%%%%%%%%% Plotting PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. For LFP 
stimPsd = squeeze(10*log10(psdData.energyValues(:,2,sfPos,oriPos,:))); 
bslPsd = squeeze(10*log10(psdData.energyValues(:,1,sfPos,oriPos,:))); 
plot(plotHand(1,1),psdData.mtFreq,mean(bslPsd/10,1),'color','k','LineWidth',1.5,'LineStyle','--');
hold(plotHand(1,1), "on");
plot(plotHand(1,1),psdData.mtFreq,mean(stimPsd/10,1),'color','k','LineWidth',1.5,'LineStyle','-');
legend(plotHand(1,1),{'baseline','stimulus'},'Box','off','FontSize',15);
xlim(plotHand(1,1), [1 150]);
set(plotHand(1,1),'NextPlot','add'); axis(plotHand(1,1),[0 150 -2 4]); set(plotHand(1,1),'TickDir','out'); box(plotHand(1,1),'off');
set(plotHand(1,1),'Xscale','log');
set(plotHand(1,1),'xTick',[10 30 70 150]); set(plotHand(1,1),'xTickLabel',[]);
set(plotHand(1,1),'yTick',[-2 0 2 4]); set(plotHand(1,1),'yTickLabel',[-2 0 2 4]);
set(plotHand(1,1),'FontSize',15);
set(plotHand(1,1),'LineWidth',1);
ylabel(plotHand(1,1), 'Raw Power (log_{10}\muV^2)','FontSize',15); 
text(60,-1.7,['\it N=' num2str(size(psdData.energyValues,1))],'Color', 'black', 'FontSize',12,'Parent',plotHand(1,1));

% Delta PSD
deltaPSD =  mean(stimPsd - bslPsd,1);
plot(deltaPlotHand(1,1),psdData.mtFreq,deltaPSD,'color','k','LineWidth',1.5,'LineStyle','-');
hold(deltaPlotHand(1,1), "on");
set(deltaPlotHand(1,1),'NextPlot','add'); 

% 2. For ECoG
modality = 'ECoG';

% Load PSD data
clear psdData
psdData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_PSD_For_All_Conditions.mat']));

% Plotting PSD
stimPsd = squeeze(10*log10(psdData.energyValues(:,2,sfPos,oriPos,:)));
bslPsd = squeeze(10*log10(psdData.energyValues(:,1,sfPos,oriPos,:)));
plot(plotHand(2,1),psdData.mtFreq,mean(bslPsd/10,1),'color',[0 0 0.7],'LineWidth',1.5,'LineStyle','--');
hold(plotHand(2,1), "on");
plot(plotHand(2,1),psdData.mtFreq,mean(stimPsd/10,1),'color',[0 0 0.7],'LineWidth',1.5,'LineStyle','-');
legend(plotHand(2,1),{'baseline','stimulus'},'Box','off','FontSize',15);
xlim(plotHand(2,1), [1 150]);
set(plotHand(2,1),'NextPlot','add'); axis(plotHand(2,1),[0 150 -2 4]); set(plotHand(2,1),'TickDir','out'); box(plotHand(2,1),'off');
set(plotHand(2,1),'Xscale','log');
set(plotHand(2,1),'xTick',[10 30 70 150]); set(plotHand(2,1),'xTickLabel',[10 30 70 150]);
set(plotHand(2,1),'yTick',[-2 0 2 4]); set(plotHand(2,1),'yTickLabel',[-2 0 2 4]);
set(plotHand(2,1),'FontSize',15);
set(plotHand(2,1),'LineWidth',1);
ylabel(plotHand(2,1), 'Raw Power (log_{10}\muV^2)','FontSize',15); xlabel(plotHand(2,1), 'Frequency (Hz)','FontSize',15); 
text(60,-1.7,['\it N=' num2str(size(psdData.energyValues,1))],'Color', 'black', 'FontSize',12,'Parent',plotHand(2,1));

% Delta PSD
deltaPSD =  mean(stimPsd - bslPsd,1);
plot(deltaPlotHand(1,1),psdData.mtFreq,deltaPSD,'color',[0 0 0.7],'LineWidth',1.5,'LineStyle','-');
hold(deltaPlotHand(1,1), "on");
xlim(deltaPlotHand(1,1), [1 150]);
set(deltaPlotHand(1,1),'NextPlot','add'); axis(deltaPlotHand(1,1),[0 150 -5 10]); set(deltaPlotHand(1,1),'TickDir','out'); box(deltaPlotHand(1,1),'off');
set(deltaPlotHand(1,1),'Xscale','log');
set(deltaPlotHand(1,1),'xTick',[10 30 70 150]); set(deltaPlotHand(1,1),'xTickLabel',[10 30 70 150]);
set(deltaPlotHand(1,1),'yTick',[-5 0 5 10]); set(deltaPlotHand(1,1),'yTickLabel',[-5 0 5 10]);
yline(deltaPlotHand(1,1),0); % plot zero line
% fast gamma
xline(deltaPlotHand(1,1),gammaRange{1}(1),'Color',[153/255 102/255 204/255],'LineWidth',1);xline(deltaPlotHand(1,1),gammaRange{1}(2),'Color',[153/255 102/255 204/255],'LineWidth',1);
% slow gamma
xline(deltaPlotHand(1,1),gammaRange{2}(1),'Color',[92/255 64/255 51/255],'LineWidth',1);xline(deltaPlotHand(1,1),gammaRange{2}(2),'Color',[92/255 64/255 51/255],'LineWidth',1);
set(deltaPlotHand(1,1),'FontSize',15);
set(deltaPlotHand(1,1),'LineWidth',1);
ylabel(deltaPlotHand(1,1), 'Change in Power (dB)','FontSize',15); xlabel(deltaPlotHand(1,1), 'Frequency (Hz)','FontSize',15);
text(2,9,'\it LFP' ,'Color', 'black', 'FontSize',12,'Parent',deltaPlotHand(1,1));
text(2,8.3,'\it ECoG' ,'Color', [0 0 0.7], 'FontSize',12,'Parent',deltaPlotHand(1,1));

annotation( 'textbox', 'String', 'Monkey 1', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight', 'Bold', 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.24,0.9,0.9,0.09], 'rotation',0 )

%% Similarly for Monkey 2

% Protocol details
monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';

%%%%%%%%%%%%%%%%%% Plotting PSD %%%%%%%%%%%%%%%%%%%%%%%%
% 1. For LFP
modality = 'LFP';

% Plot Handles
psdPlotsM2 = getPlotHandles(2,1,[0.53 0.34 0.15 0.6],0.01,0.03);
deltaPsdPlotsM2 = getPlotHandles(1,1,[0.74 0.34 0.15 0.6],0.01,0.01);
plotHand = psdPlotsM2;
deltaPlotHand = deltaPsdPlotsM2;

% load PSD data
psdData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_spike' modality '_PSD_For_All_Conditions.mat']));

%%%%%%%%%%%%%%%%%%%%% Plotting PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. For LFP 
stimPsd = squeeze(10*log10(psdData.energyValues(:,2,sfPos,oriPos,:)));
bslPsd = squeeze(10*log10(psdData.energyValues(:,1,sfPos,oriPos,:)));
plot(plotHand(1,1),psdData.mtFreq,mean(bslPsd/10,1),'color','k','LineWidth',1.5,'LineStyle','--');
hold(plotHand(1,1), "on");
plot(plotHand(1,1),psdData.mtFreq,mean(stimPsd/10,1),'color','k','LineWidth',1.5,'LineStyle','-');
% legend(plotHand(1,1),{'baseline','stimulus'},'Box','off','FontSize',15);
xlim(plotHand(1,1), [1 150]);
set(plotHand(1,1),'NextPlot','add'); axis(plotHand(1,1),[0 150 -2 4]); set(plotHand(1,1),'TickDir','out'); box(plotHand(1,1),'off');
set(plotHand(1,1),'Xscale','log');
set(plotHand(1,1),'xTick',[10 30 70 150]); set(plotHand(1,1),'xTickLabel',[]);
set(plotHand(1,1),'yTick',[-2 0 2 4]); set(plotHand(1,1),'yTickLabel',[-2 0 2 4]);
set(plotHand(1,1),'FontSize',15);
set(plotHand(1,1),'LineWidth',1);
ylabel(plotHand(1,1), 'Raw Power (log_{10}\muV^2)','FontSize',15); 
text(60,-1.7,['\it N=' num2str(size(psdData.energyValues,1))],'Color', 'black', 'FontSize',12,'Parent',plotHand(1,1));
text(20,3,'\it LFP' ,'Color', 'black', 'FontSize',12,'Parent',plotHand(1,1));

% Delta PSD
deltaPSD =  mean(stimPsd - bslPsd,1);
plot(deltaPlotHand(1,1),psdData.mtFreq,deltaPSD,'color','k','LineWidth',1.5,'LineStyle','-');
hold(deltaPlotHand(1,1), "on");
set(deltaPlotHand(1,1),'NextPlot','add'); 

% 2. For ECoG
modality = 'ECoG';

% load PSD data
clear psdData
psdData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_PSD_For_All_Conditions.mat']));

% Plotting PSD
stimPsd = squeeze(10*log10(psdData.energyValues(:,2,sfPos,oriPos,:)));
bslPsd = squeeze(10*log10(psdData.energyValues(:,1,sfPos,oriPos,:)));
plot(plotHand(2,1),psdData.mtFreq,mean(bslPsd/10,1),'color',[0 0 0.7],'LineWidth',1.5,'LineStyle','--');
hold(plotHand(2,1), "on");
plot(plotHand(2,1),psdData.mtFreq,mean(stimPsd/10,1),'color',[0 0 0.7],'LineWidth',1.5,'LineStyle','-');
xlim(plotHand(2,1), [1 150]);
set(plotHand(2,1),'NextPlot','add'); axis(plotHand(2,1),[0 150 -2 4]); set(plotHand(2,1),'TickDir','out'); box(plotHand(2,1),'off');
set(plotHand(2,1),'Xscale','log');
set(plotHand(2,1),'xTick',[10 30 70 150]); set(plotHand(2,1),'xTickLabel',[10 30 70 150]);
set(plotHand(2,1),'yTick',[-2 0 2 4]); set(plotHand(2,1),'yTickLabel',[-2 0 2 4]);
set(plotHand(2,1),'FontSize',15);
set(plotHand(2,1),'LineWidth',1);
ylabel(plotHand(2,1), 'Raw Power (log_{10}\muV^2)','FontSize',15); xlabel(plotHand(2,1), 'Frequency (Hz)','FontSize',15); 
text(60,-1.7,['\it N=' num2str(size(psdData.energyValues,1))],'Color', 'black', 'FontSize',12,'Parent',plotHand(2,1));
text(20,3,'\it ECoG' ,'Color', [0.5 0.5 0.5], 'FontSize',12,'Parent',plotHand(2,1));

% Delta PSD
deltaPSD =  mean(stimPsd - bslPsd,1);
plot(deltaPlotHand(1,1),psdData.mtFreq,deltaPSD,'color',[0 0 0.7],'LineWidth',1.5,'LineStyle','-');
hold(deltaPlotHand(1,1), "on");
xlim(deltaPlotHand(1,1), [1 150]);
set(deltaPlotHand(1,1),'NextPlot','add'); axis(deltaPlotHand(1,1),[0 150 -5 15]); set(deltaPlotHand(1,1),'TickDir','out'); box(deltaPlotHand(1,1),'off');
set(deltaPlotHand(1,1),'Xscale','log');
set(deltaPlotHand(1,1),'xTick',[10 30 70 150]); set(deltaPlotHand(1,1),'xTickLabel',[10 30 70 150]);
set(deltaPlotHand(1,1),'yTick',[-5 0 5 10 15]); set(deltaPlotHand(1,1),'yTickLabel',[-5 0 5 10 15]);
% plot zero line
yline(deltaPlotHand(1,1),0)
% fast gamma
xline(deltaPlotHand(1,1),gammaRange{1}(1),'Color',[153/255 102/255 204/255],'LineWidth',1);xline(deltaPlotHand(1,1),gammaRange{1}(2),'Color',[153/255 102/255 204/255],'LineWidth',1);
% slow gamma
xline(deltaPlotHand(1,1),gammaRange{2}(1),'Color',[92/255 64/255 51/255],'LineWidth',1);xline(deltaPlotHand(1,1),gammaRange{2}(2),'Color',[92/255 64/255 51/255],'LineWidth',1);
set(deltaPlotHand(1,1),'FontSize',15);
set(deltaPlotHand(1,1),'LineWidth',1);
ylabel(deltaPlotHand(1,1), 'Change in Power (dB)','FontSize',15); xlabel(deltaPlotHand(1,1), 'Frequency (Hz)','FontSize',15);
text(2,14,'\it LFP' ,'Color', 'black', 'FontSize',12,'Parent',deltaPlotHand(1,1));
text(2,13,'\it ECoG' ,'Color', [0 0 0.7], 'FontSize',12,'Parent',deltaPlotHand(1,1));

annotation( 'textbox', 'String', 'Monkey 2', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight', 'Bold', 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.675,0.9,0.9,0.09], 'rotation',0 )



