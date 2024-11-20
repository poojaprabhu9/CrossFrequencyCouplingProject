% Figure 4A-F: Plotting for Phase-Amplitude Coupling (PAC) for SF='all' and Ori=0:22.5:180 conditions
% for Local Field Potential (LFP), under the case, D0 after MP (single electrode LFP after Matching Pursuit(MP))

% The line curves are generated for two different frequency bands (Slow Gamma, and Fast Gamma)

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms

clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Electrode choices
modality = 'LFP'; 

% Datapath
folderSave = 'savedData';

% colormap for PAC Plots
colormap jet

% SF = 'all' and Ori all conditions
sfPos = 6; oriPos = [1 2 3 4 5 6 7 8];

% gamma range
gammaRange = {[40 70], [25 40]}; % Fast Gamma and Slow Gamma

% Frequency range for three bounding boxes-Fast gamma and slow gamma, respectively
stimFreqWin = {[gammaRange{1}(1) gammaRange{1}(2) 167 487], [gammaRange{2}(1) gammaRange{2}(2) 87 157]};

% Color for two boxes-Fast gamma and slow gamma, respectively
stimFreqWinColors = [153/255 102/255 204/255;  92/255 64/255 51/255];

% For better visualization: splitting PAC plot tto two parts: 2 to 150 Hz and 150Hz to 500Hz
partitionFreq = [7 157 157 487]; 

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];
pacMethod = 'klmi'; 
filterName = 'fir'; 
nSurrogates = 0; % Number of surrogates
useMPFlag = 1; % Particular to MP filtered data
sVarName = 'sf';
electrodeDistanceVal = '0';
titleNames = {'0^{o}', '22.5^{o}', '45^{o}', '67.5^{o}', '90^{o}', '112.5^{o}', '135^{o}', '157.5^{o}'};

%% Monkey 1

% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol

% Plot Handles
pacPlotsM1 = getPlotHandles(2,8,[0.05 0.38 0.40 0.53],0.005,0);
linePlotsM1 = getPlotHandles(1,2,[0.05 0.09 0.40 0.18],0.05,0.01);
pacPlotHand = pacPlotsM1;
linePlotHand = linePlotsM1;

% Color bar limits for 150-500 Hz and 2-150 Hz PAC plots
climVals1 = 0.001; climVals2 = 0.0005;

% load PAC
for i = 1:length(oriPos)
    clear tmpData
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos(i)) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
    freqIdx = find(tmpData.centerAmpFreq == partitionFreq(2));

    % 150 to 500Hz
    pcolor(pacPlotHand(1,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(freqIdx:end),squeeze(mean(tmpData.pac(:,2,freqIdx:end,:),1))); 
    shading(pacPlotHand(1,i),'interp'); box(pacPlotHand(1,i),'off');
    hold(pacPlotHand(1,i), "on");
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,1}(1) stimFreqWin{1,1}(2)], [stimFreqWin{1,1}(3) stimFreqWin{1,1}(4)], stimFreqWinColors(1,:), 2, '-', 'B');
    colormap(pacPlotHand(1,i),'jet');
    clim(pacPlotHand(1,i), [0 climVals1]);
    set(pacPlotHand(1,i),'Xscale','log');
    set(pacPlotHand(1,i),'Yscale','log');
    set(pacPlotHand(1,i),'FontSize',10);

    % 2 to 150Hz
    pcolor(pacPlotHand(2,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(1:freqIdx),squeeze(mean(tmpData.pac(:,2,1:freqIdx,:),1))); 
    shading(pacPlotHand(2,i),'interp'); box(pacPlotHand(2,i),'off');
    hold(pacPlotHand(2,i), "on");
    makeBox(pacPlotHand(2,i), [stimFreqWin{1,2}(1) stimFreqWin{1,2}(2)], [stimFreqWin{1,2}(3) stimFreqWin{1,2}(4)], stimFreqWinColors(2,:), 2, '-', 'B');
    colormap(pacPlotHand(2,i),'jet');
    clim(pacPlotHand(2,i), [0 climVals2]);
    set(pacPlotHand(2,i),'Xscale','log');
    set(pacPlotHand(2,i),'Yscale','log');
    set(pacPlotHand(2,i),'FontSize',10);
    set(pacPlotHand(2,i),'xTick',[10 30 70]); set(pacPlotHand(2,i),'xTickLabel',[10 30 70]); 
    title(pacPlotHand(1,i), titleNames{i}, 'FontSize', 12);

    if i==1
        set(pacPlotHand(1,i),'yTick',[150 250 350 450]); set(pacPlotHand(1,i),'yTickLabel',[ 150 250 350 450]); 
        set(pacPlotHand(2,i),'yTick',[ 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[  50 100 150]); 

        xlabel(pacPlotHand(2,i),'Phase Frequency (Hz)','Fontsize',15);
        ylabel(pacPlotHand(2,i),'Amplitude Frequency (Hz)','Fontsize',15);
    else
        set(pacPlotHand(1,i),'yTick',[150 250 350 450]); set(pacPlotHand(1,i),'yTickLabel',[]); 
        set(pacPlotHand(2,i),'yTick',[ 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[]); 
    end
end

cb=colorbar(pacPlotHand(1,8));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.028 0 0 -0.18];
cb.Ticks = [0 climVals1];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',12,'Rotation',270);

cb=colorbar(pacPlotHand(2,8));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.028 0 0 -0.18];
cb.Ticks = [0 climVals2];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',12,'Rotation',270);

%%%%%%%% Plotting tuning curve for PAC values for each gamma range %%%%%%%%%%%%%
% initialise
meanVals = cell(1, length(stimFreqWin));
semVals = cell(1, length(stimFreqWin));
xtickLabels = [0 22.5 45 67.5 90 112.5 135 157.5];

for nWin = 1:length(stimFreqWin)
    if mod(stimFreqWin{1,nWin}(1),4)==1
        xFreq1=stimFreqWin{1,nWin}(1)-1;
    elseif mod(stimFreqWin{1,nWin}(1),4)==2
        xFreq1=stimFreqWin{1,nWin}(1)+2;
    elseif mod(stimFreqWin{1,nWin}(1),4)==3
        xFreq1=stimFreqWin{1,nWin}(1)+1;
    elseif mod(stimFreqWin{1,nWin}(1),4)==0
        xFreq1=stimFreqWin{1,nWin}(1);
    end

    if mod(stimFreqWin{1,nWin}(2),4)==1
        xFreq2=stimFreqWin{1,nWin}(2)-1;
    elseif mod(stimFreqWin{1,nWin}(2),4)==2
        xFreq2=stimFreqWin{1,nWin}(2)+2;
    elseif mod(stimFreqWin{1,nWin}(2),4)==3
        xFreq2=stimFreqWin{1,nWin}(2)+1;
    elseif mod(stimFreqWin{1,nWin}(2),4)==0
        xFreq2=stimFreqWin{1,nWin}(2);
    end
    xFreqPos = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));
    yFreqPos = intersect(find(tmpData.centerAmpFreq>= stimFreqWin{1,nWin}(3)),find(tmpData.centerAmpFreq<= stimFreqWin{1,nWin}(4)));

    for i = 1:length(oriPos)
        clear tmpData
        tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos(i)) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
        avgPacElec = squeeze(mean(tmpData.pac(:,2,yFreqPos,xFreqPos),[3,4]));
        meanVals{1,nWin}(i,1) = mean(avgPacElec,'all');
        pacPerElec = squeeze(tmpData.pac(:,2,yFreqPos,xFreqPos));
        semVals{1,nWin}(i,1) = getSEMedian(pacPerElec(:),1000);
    end
    errorbar(linePlotHand(1,1),xtickLabels,meanVals{1,nWin}(:,1),semVals{1,nWin}(:,1),'color',stimFreqWinColors(nWin,:),'lineWidth',1.5,'Marker','o','CapSize',0,'MarkerSize',3,'MarkerFaceColor','black','MarkerEdgeColor','black');
    set(linePlotHand(1,1),'NextPlot','add');  
end

set(linePlotHand(1,1),'TickDir','out'); 
set(linePlotHand(1,1),'xTick',xtickLabels); set(linePlotHand(1,1),'xTickLabel',xtickLabels);
box(linePlotHand(1,1),'off');
set(linePlotHand(1,1),'FontSize',10);
set(linePlotHand(1,1),'LineWidth',1);
xlabel(linePlotHand(1,1),'Orientation (degrees)','Fontsize',15);
ylabel(linePlotHand(1,1),'Averaged KLMI','Fontsize',15);
set(linePlotHand(1,1),'xTick',[0 45 90 135]); set(linePlotHand(1,1),'xTickLabel',[0 45 90 135]);
text(120,7.6e-4,['\it N=' num2str(size(tmpData.pac,1))],'Color', 'black', 'Fontsize',12,'Parent',linePlotHand(1,1));

%%%%%%%% Plotting tuning curve for PSD values for each gamma range %%%%%%%%%%%%%
% Signal processing choices
tapers = [2 3];
psdData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_spike' modality '_PSD_For_All_Conditions.mat']));
stimPsd = squeeze(10*log10(psdData.energyValues(:,2,sfPos,oriPos,:)));
bslPsd = squeeze(10*log10(psdData.energyValues(:,1,sfPos,oriPos,:)));

meanPsdVals = cell(1, length(stimFreqWin));
semPsdVals = cell(1, length(stimFreqWin));

for nWin = 1:length(stimFreqWin)
    gammaPos = intersect(find(psdData.mtFreq>=gammaRange{1,nWin}(1)),find(psdData.mtFreq<=gammaRange{1,nWin}(2)));
    for i = 1:length(oriPos)
        stimPsdPerCond = squeeze(stimPsd(:,i,:)); 
        bslPsdPerCond = squeeze(bslPsd(:,i,:));

        % Delta PSD
        deltaPSD = stimPsdPerCond(:,gammaPos) - bslPsdPerCond(:,gammaPos);
        avgPsdElec = mean(deltaPSD,2);
        meanPsdVals{1,nWin}(i,1) = mean(avgPsdElec,'all');
        semPsdVals{1,nWin}(i,1) = getSEMedian(deltaPSD(:),1000);
    end
    errorbar(linePlotHand(1,2),xtickLabels,meanPsdVals{1,nWin}(:,1),semPsdVals{1,nWin}(:,1),'color',stimFreqWinColors(nWin,:),'lineWidth',1.5,'Marker','o','CapSize',0,'MarkerSize',3,'MarkerFaceColor','black','MarkerEdgeColor','black');
    set(linePlotHand(1,2),'NextPlot','add');
end

set(linePlotHand(1,2),'TickDir','out'); 
set(linePlotHand(1,2),'xTick',xtickLabels); set(linePlotHand(1,2),'xTickLabel',xtickLabels);
box(linePlotHand(1,2),'off');
set(linePlotHand(1,2),'FontSize',10);
set(linePlotHand(1,2),'LineWidth',1);
xlabel(linePlotHand(1,2),'Orientation (degrees)','Fontsize',15);
ylabel(linePlotHand(1,2),'Averaged Power (dB)','Fontsize',15);

ylim(linePlotHand(1,2),[3 9])
set(linePlotHand(1,2),'yTick',[6 9]); set(linePlotHand(1,2),'yTickLabel',[6 9]);
set(linePlotHand(1,2),'xTick',[0 45 90 135]); set(linePlotHand(1,2),'xTickLabel',[0 45 90 135]);
text(120,8.6,['\it N=' num2str(size(tmpData.pac,1))],'Color', 'black', 'Fontsize',12,'Parent',linePlotHand(1,2));

annotation( 'textbox', 'String', 'Monkey 1', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight', 'Bold','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.24,0.3,0.9,0.7] )

%%%%%%%%% Correlation %%%%%%%%%
[R,P] = corr(meanVals{1, 1},meanPsdVals{1, 1},'type','Spearman');
disp(['correlation between FG_pac and FG_pow=' num2str(R) 'with p-val=' num2str(P)]);

[R,P] = corr(meanVals{1, 1},meanPsdVals{1, 2},'type','Spearman');
disp(['correlation between FG_pac and SG_pow=' num2str(R) 'with p-val=' num2str(P)]);

[R,P] = corr(meanVals{1, 2},meanPsdVals{1, 2},'type','Spearman');
disp(['correlation between SG_pac and SG_pow=' num2str(R) 'with p-val=' num2str(P)]);

[R,P] = corr(meanVals{1, 2},meanPsdVals{1, 1},'type','Spearman');
disp(['correlation between SG_pac and FG_pow=' num2str(R) 'with p-val=' num2str(P)]);

%% Monkey 2
monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001';% SF-Ori protocol

% Plot Handles
pacPlotsM2 = getPlotHandles(2,8,[0.55 0.38 0.40 0.53],0.005,0);
linePlotsM2 = getPlotHandles(1,2,[0.55 0.09 0.40 0.18],0.05,0.01);
pacPlotHand = pacPlotsM2;
linePlotHand = linePlotsM2;

% Color bar limits for 150-500 Hz and 2-150 Hz PAC plots
climVals1 = 0.0005; climVals2 = 0.00025;

% load PAC
for i = 1:length(oriPos)
    clear tmpData
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos(i)) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
    freqIdx = find(tmpData.centerAmpFreq == partitionFreq(2)); 

    % 150 to 500Hz
    pcolor(pacPlotHand(1,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(freqIdx:end),squeeze(mean(tmpData.pac(:,2,freqIdx:end,:),1))); 
    shading(pacPlotHand(1,i),'interp'); box(pacPlotHand(1,i),'off');
    hold(pacPlotHand(1,i), "on");
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,1}(1) stimFreqWin{1,1}(2)], [stimFreqWin{1,1}(3) stimFreqWin{1,1}(4)], stimFreqWinColors(1,:), 2, '-', 'B');
    colormap(pacPlotHand(1,i),'jet');
    clim(pacPlotHand(1,i), [0 climVals1]);
    set(pacPlotHand(1,i),'Xscale','log');
    set(pacPlotHand(1,i),'Yscale','log');
    set(pacPlotHand(1,i), 'FontSize',10);

    % 2 to 150Hz
    pcolor(pacPlotHand(2,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(1:freqIdx),squeeze(mean(tmpData.pac(:,2,1:freqIdx,:),1))); 
    shading(pacPlotHand(2,i),'interp'); box(pacPlotHand(2,i),'off');
    hold(pacPlotHand(2,i), "on");
    makeBox(pacPlotHand(2,i), [stimFreqWin{1,2}(1) stimFreqWin{1,2}(2)], [stimFreqWin{1,2}(3) stimFreqWin{1,2}(4)], stimFreqWinColors(2,:), 2, '-', 'B');
    colormap(pacPlotHand(2,i),'jet');
    clim(pacPlotHand(2,i), [0 climVals2]);
    set(pacPlotHand(2,i),'Xscale','log');
    set(pacPlotHand(2,i),'Yscale','log');
    set(pacPlotHand(2,i), 'FontSize',10);
    set(pacPlotHand(2,i),'xTick',[10 30 70]); set(pacPlotHand(2,i),'xTickLabel',[10 30 70]); 
    title(pacPlotHand(1,i), titleNames{i}, 'FontSize', 12);

    if i==1
        set(pacPlotHand(1,i),'yTick',[150 250 350 450]); set(pacPlotHand(1,i),'yTickLabel',[ 150 250 350 450]); 
        set(pacPlotHand(2,i),'yTick',[ 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[ 50 100 150]); 

        xlabel(pacPlotHand(2,i),'Phase Frequency (Hz)','Fontsize',15);
        ylabel(pacPlotHand(2,i),'Amplitude Frequency (Hz)','Fontsize',15);
    else
        set(pacPlotHand(1,i),'yTick',[150 250 350 450]); set(pacPlotHand(1,i),'yTickLabel',[]); 
        set(pacPlotHand(2,i),'yTick',[ 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[]); 
    end

end
cb=colorbar(pacPlotHand(1,8));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.028 0 0 -0.18];
cb.Ticks = [0 climVals1];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',12,'Rotation',270);

cb=colorbar(pacPlotHand(2,8));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.028 0 0 -0.18];
cb.Ticks = [0 climVals2];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',12,'Rotation',270);

%%%%%%%%%%%% Plotting tuning curve for PAC values for each gamma range %%%%%%%%%%%%%%%%%%%%%%%
% initialise
meanVals = cell(1, length(stimFreqWin));
semVals = cell(1, length(stimFreqWin));
xtickLabels = [0 22.5 45 67.5 90 112.5 135 157.5];

for nWin = 1:length(stimFreqWin)
    if mod(stimFreqWin{1,nWin}(1),4)==1
        xFreq1=stimFreqWin{1,nWin}(1)-1;
    elseif mod(stimFreqWin{1,nWin}(1),4)==2
        xFreq1=stimFreqWin{1,nWin}(1)+2;
    elseif mod(stimFreqWin{1,nWin}(1),4)==3
        xFreq1=stimFreqWin{1,nWin}(1)+1;
    elseif mod(stimFreqWin{1,nWin}(1),4)==0
        xFreq1=stimFreqWin{1,nWin}(1);
    end

    if mod(stimFreqWin{1,nWin}(2),4)==1
        xFreq2=stimFreqWin{1,nWin}(2)-1;
    elseif mod(stimFreqWin{1,nWin}(2),4)==2
        xFreq2=stimFreqWin{1,nWin}(2)+2;
    elseif mod(stimFreqWin{1,nWin}(2),4)==3
        xFreq2=stimFreqWin{1,nWin}(2)+1;
    elseif mod(stimFreqWin{1,nWin}(2),4)==0
        xFreq2=stimFreqWin{1,nWin}(2);
    end
    xFreqPos = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));
    yFreqPos = intersect(find(tmpData.centerAmpFreq>= stimFreqWin{1,nWin}(3)),find(tmpData.centerAmpFreq<= stimFreqWin{1,nWin}(4)));

    for i = 1:length(oriPos)
        clear tmpData
        tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sfPos) '_o' num2str(oriPos(i)) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
        avgPacElec = squeeze(mean(tmpData.pac(:,2,yFreqPos,xFreqPos),[3,4]));
        meanVals{1,nWin}(i,1) = mean(avgPacElec,'all');
        pacPerElec = squeeze(tmpData.pac(:,2,yFreqPos,xFreqPos));
        semVals{1,nWin}(i,1) = getSEMedian(pacPerElec(:),1000);
    end
    errorbar(linePlotHand(1,1),xtickLabels,meanVals{1,nWin}(:,1),semVals{1,nWin}(:,1),'color',stimFreqWinColors(nWin,:),'lineWidth',1.5,'Marker','o','CapSize',0,'MarkerSize',3,'MarkerFaceColor','black','MarkerEdgeColor','black');
    set(linePlotHand(1,1),'NextPlot','add');
end

set(linePlotHand(1,1),'TickDir','out'); 
set(linePlotHand(1,1),'xTick',xtickLabels); set(linePlotHand(1,1),'xTickLabel',xtickLabels);
box(linePlotHand(1,1),'off');
set(linePlotHand(1,1),'FontSize',10);
set(linePlotHand(1,1),'LineWidth',1);
xlabel(linePlotHand(1,1),'Orientation (degrees)','Fontsize',15);
ylabel(linePlotHand(1,1),'Averaged KLMI','Fontsize',15);

set(linePlotHand(1,1),'xTick',[0 45 90 135]); set(linePlotHand(1,1),'xTickLabel',[0 45 90 135]);
text(120,5.6e-4,['\it N=' num2str(size(tmpData.pac,1))],'Color', 'black', 'Fontsize',12,'Parent',linePlotHand(1,1));

%%%%%%%% Plotting tuning curve for PSD values for each gamma range %%%%%%%%%%%%%
% Signal processing choices
tapers = [2 3];
clear psdData
psdData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_spike' modality '_PSD_For_All_Conditions.mat']));
stimPsd = squeeze(10*log10(psdData.energyValues(:,2,sfPos,oriPos,:)));
bslPsd = squeeze(10*log10(psdData.energyValues(:,1,sfPos,oriPos,:)));

meanPsdVals = cell(1, length(stimFreqWin));
semPsdVals = cell(1, length(stimFreqWin));

for nWin = 1:length(stimFreqWin)
    gammaPos = intersect(find(psdData.mtFreq>=gammaRange{1,nWin}(1)),find(psdData.mtFreq<=gammaRange{1,nWin}(2)));
    for i = 1:length(oriPos)
        stimPsdPerCond = squeeze(stimPsd(:,i,:)); 
        bslPsdPerCond = squeeze(bslPsd(:,i,:));

        % Delta PSD
        deltaPSD = stimPsdPerCond(:,gammaPos) - bslPsdPerCond(:,gammaPos);
        avgPsdElec = mean(deltaPSD,2);
        meanPsdVals{1,nWin}(i,1) = mean(avgPsdElec,'all');
        semPsdVals{1,nWin}(i,1) = getSEMedian(deltaPSD(:),1000);
    end
    errorbar(linePlotHand(1,2),xtickLabels,meanPsdVals{1,nWin}(:,1),semPsdVals{1,nWin}(:,1),'color',stimFreqWinColors(nWin,:),'lineWidth',1.5,'Marker','o','CapSize',0,'MarkerSize',3,'MarkerFaceColor','black','MarkerEdgeColor','black');
    set(linePlotHand(1,2),'NextPlot','add');
end

set(linePlotHand(1,2),'TickDir','out'); 
set(linePlotHand(1,2),'xTick',xtickLabels); set(linePlotHand(1,2),'xTickLabel',xtickLabels);
box(linePlotHand(1,2),'off');
set(linePlotHand(1,2),'FontSize',10);
set(linePlotHand(1,2),'LineWidth',1);
xlabel(linePlotHand(1,2),'Orientation (degrees)','Fontsize',15);
ylabel(linePlotHand(1,2),'Averaged Power (dB)','Fontsize',15);

ylim(linePlotHand(1,2),[4 8])
set(linePlotHand(1,2),'yTick',[4 6 8]); set(linePlotHand(1,2),'yTickLabel',[4 6 8]);
set(linePlotHand(1,2),'xTick',[0 45 90 135]); set(linePlotHand(1,2),'xTickLabel',[0 45 90 135]);
text(120,7.7,['\it N=' num2str(size(tmpData.pac,1))],'Color', 'black', 'Fontsize',12,'Parent',linePlotHand(1,2));

annotation( 'textbox', 'String', 'Monkey 2', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight', 'Bold','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.74,0.3,0.9,0.7] )

%%%%%%%%% Correlation Values %%%%%%%%%
[R,P] = corr(meanVals{1, 1},meanPsdVals{1, 1},'type','Spearman');
disp(['correlation between FG_pac and FG_pow=' num2str(R) 'with p-val=' num2str(P)]);

[R,P] = corr(meanVals{1, 1},meanPsdVals{1, 2},'type','Spearman');
disp(['correlation between FG_pac and SG_pow=' num2str(R) 'with p-val=' num2str(P)]);

[R,P] = corr(meanVals{1, 2},meanPsdVals{1, 2},'type','Spearman');
disp(['correlation between SG_pac and SG_pow=' num2str(R) 'with p-val=' num2str(P)]);

[R,P] = corr(meanVals{1, 2},meanPsdVals{1, 1},'type','Spearman');
disp(['correlation between SG_pac and FG_pow=' num2str(R) 'with p-val=' num2str(P)]);
