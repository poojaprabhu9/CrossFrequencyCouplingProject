% Figure 3A-D: Plotting for Phase-Amplitude Coupling (PAC) for 'all-all'
% conditions for Local Field Potential (LFP), under three cases,
% 1. D0 (single electrode LFP)
% 2. D0 after MP (single electrode LFP after Matching Pursuit(MP))
% 3. D400 (inter-electrode within 400 micrometer distance)
% For all three cases the violin plots are generated for three different
% frequency bands (Low Frequency, Slow Gamma, and Fast Gamma)

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms
% gramm data visualization toolbox: https://github.com/piermorel/gramm

clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Electrode choices
modality = 'LFP'; 

% Data path
folderSave = 'savedData';

% Colormap for PAC Plots
colormap jet

% SF and Ori position for 'all-all' conditions
sfPos = 6; oriPos = 9;

% Fast Gamma and Slow Gamma range, respectively
gammaRange = {[40 70], [25 40]}; % Fast Gamma and Slow Gamma; Chosen as per the peaks visible in PSD in both the Monkeys

% Frequency range for three bounding boxes-Fast gamma, slow gamma and low
% frequency, respectively
stimFreqWin = {[gammaRange{1}(1) gammaRange{1}(2) 167 487], [gammaRange{2}(1) gammaRange{2}(2) 87 157], [4 8 17 67] };

% For better visualization: splitting PAC plot tto two parts: 2 to 150 Hz and 150Hz to 500Hz
partitionFreq = [7 157 157 487];

% Color for three boxes-Fast gamma, slow gamma and low
% frequency, respectively
stimFreqWinColors = [153/255 102/255 204/255;  92/255 64/255 51/255; 0 0 0];

% labels for each frequency windows
labelNames = { {'Fast Gamma','with 150-500 Hz'}, {'Slow Gamma','with 80-150 Hz'}, {'Low Frequency(<10 Hz)','with 20-70 Hz'}};

% labels needed to create the structure
xLabels = {'D0', 'D0_afterMP', 'D_400'};

% Color vals for three cases (D0, D0 after MP and D400)
colorVals = [1 0 1; 0.2 0.8 0.5; 1 0.6824 0.3176]; 

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];
pacMethod = 'klmi'; 
filterName = 'fir'; 
nSurrogates = 0; % Number of surrogates
useMPFlag = {'0','1','0'}; % Particular to MP filtered data
sVarName = 'sf';
electrodeDistanceVal = {'0','0','400'};
titleNames = {'D_{0}', 'D_{0} after MP', 'D_{400}'};

%% Monkey 1

% Protocol details
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol

% Plot Handles
pacPlotsM1 = getPlotHandles(2,3,[0.05 0.4 0.40 0.53],0.005,0);
pacPlotHand = pacPlotsM1;

% Color bar limits for 150-500 Hz and 2-150 Hz PAC plots
climVals1 = 0.001; climVals2 = 0.0005;

% load PAC
for i = 1:length(electrodeDistanceVal)
    clear tmpData
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal{i} '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag{i}) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
    freqIdx = find(tmpData.centerAmpFreq == partitionFreq(2));
    %150 to 500Hz
    pcolor(pacPlotHand(1,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(freqIdx:end),squeeze(mean(tmpData.pac(:,2,freqIdx:end,:),1))); 
    shading(pacPlotHand(1,i),'interp'); box(pacPlotHand(1,i),'off');
    hold(pacPlotHand(1,i), "on");
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,1}(1) stimFreqWin{1,1}(2)], [stimFreqWin{1,1}(3) stimFreqWin{1,1}(4)], stimFreqWinColors(1,:), 2, '-', 'B');
    colormap(pacPlotHand(1,i),'jet');
    clim(pacPlotHand(1,i), [0 climVals1]);
    set(pacPlotHand(1,i),'Xscale','log');
    set(pacPlotHand(1,i),'Yscale','log');
    set(pacPlotHand(1,i),'FontSize',12);

    % 2 to 150Hz
    pcolor(pacPlotHand(2,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(1:freqIdx),squeeze(mean(tmpData.pac(:,2,1:freqIdx,:),1))); 
    shading(pacPlotHand(2,i),'interp'); box(pacPlotHand(2,i),'off');
    hold(pacPlotHand(2,i), "on");
    makeBox(pacPlotHand(2,i), [stimFreqWin{1,2}(1) stimFreqWin{1,2}(2)], [stimFreqWin{1,2}(3) stimFreqWin{1,2}(4)], stimFreqWinColors(2,:), 2, '-', 'B');
    makeBox(pacPlotHand(2,i), [stimFreqWin{1,3}(1) stimFreqWin{1,3}(2)], [stimFreqWin{1,3}(3) stimFreqWin{1,3}(4)], stimFreqWinColors(3,:), 2, '-', 'B');
    colormap(pacPlotHand(2,i),'jet');
    clim(pacPlotHand(2,i), [0 climVals2]);
    set(pacPlotHand(2,i),'Xscale','log');
    set(pacPlotHand(2,i),'Yscale','log');
    set(pacPlotHand(2,i),'FontSize',12);
    set(pacPlotHand(2,i),'xTick',[10 30 70]); set(pacPlotHand(2,i),'xTickLabel',[10 30 70]); 
    title(pacPlotHand(1,i), [titleNames{i} '(N=' num2str(size(tmpData.pac,1)) ')'], 'Color', colorVals(i,:), 'FontSize', 12);

    if i==1
        set(pacPlotHand(1,i),'yTick',[150 250 350 450]); set(pacPlotHand(1,i),'yTickLabel',[ 150 250 350 450]); 
        set(pacPlotHand(2,i),'yTick',[ 10 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[ 10 50 100 150]); 
        xlabel(pacPlotHand(2,i),'Phase Frequency (Hz)','Fontsize',12);
        ylabel(pacPlotHand(2,i),'Amplitude Frequency (Hz)','Fontsize',12);
    else
        set(pacPlotHand(1,i),'yTick',[150 250 350 450]); set(pacPlotHand(1,i),'yTickLabel',[]); 
        set(pacPlotHand(2,i),'yTick',[ 10 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[]); 
    end
end

cb=colorbar(pacPlotHand(1,3));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.065 0 0 -0.18];
cb.Ticks = [0 climVals1];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',10,'Rotation',270);

cb=colorbar(pacPlotHand(2,3));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.0705 0 0 -0.18];
cb.Ticks = [0 climVals2];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',10,'Rotation',270);

set(pacPlotHand,'TickLength',[0.03 0.03]);
set(pacPlotHand(2,1),'xTick',[ 5 10 30 70]); set(pacPlotHand(2,1),'xTickLabel',[ 5 10 30 70]);  
set(pacPlotHand(2,2),'xTick',[ 5 10 30 70]); set(pacPlotHand(2,2),'xTickLabel',[ 5 10 30 70]); 
set(pacPlotHand(2,3),'xTick',[ 5 10 30 70]); set(pacPlotHand(2,3),'xTickLabel',[ 5 10 30 70]); 

set(pacPlotHand(1,1),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,1),'yTickLabel',[5 10 20 30 40 70 100 300 400]); 
set(pacPlotHand(1,2),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,2),'yTickLabel',[]);
set(pacPlotHand(1,3),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,3),'yTickLabel',[]);
xlabel(pacPlotHand(1,1),'Phase Frequency (Hz)','Fontsize',12);
ylabel(pacPlotHand(1,1),'Amplitude Frequency (Hz)','Fontsize',12);
annotation( 'textbox', 'String', 'Monkey 1', 'Color', 'black', ...
            'FontSize', 14,'FontWeight','Bold','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.22,0.3,0.9,0.7] )

% total number of frequency windows
numFreqWin = length(stimFreqWin) ;

% initialise
pacVals = cell(1, numFreqWin);
pacLabelsPerWin = cell(1, numFreqWin);
pacValsPerDistPerWin = cell(1, numFreqWin);

% axis position
axisPos = {[0.30 0.02 0.17 0.30], [0.16 0.02 0.17 0.30], [0.01 0.02 0.18 0.30]};

for nWin = 1:numFreqWin
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

    pacElec = [];   pacElecTmp = []; pacLabels = []; clear pacLabelsTmp; clear pacPerElecDist;
    for i = 1:length(electrodeDistanceVal)
        tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal{i} '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag{i}) '.mat']), 'pac', 'centerPhaseFreq', 'centerAmpFreq');
        pacPerElec = squeeze(mean(tmpData.pac(:,2,yFreqPos,xFreqPos),[3,4]));
        % PAC for each box
        pacPerElecDist(i,:,:) = squeeze(mean(tmpData.pac(:,2,yFreqPos,xFreqPos),1));
        pacElecTmp = pacPerElec(:);
        for nCh = 1:length(pacElecTmp)
            pacLabelsTmp{nCh,1} = xLabels{i};
        end
        pacLabels = cat(1,pacLabels,pacLabelsTmp);
        pacElec = cat(1,pacElec,pacElecTmp);       
    end
    pacVals{1,nWin} = pacElec;
    pacLabelsPerWin{1,nWin} = pacLabels;
    pacValsPerDistPerWin{1,nWin} = pacPerElecDist;
    
end

%%%%%% Stats and Reduction Percentage %%%%%%%%%%%%%%%%%
for nWin = 1:numFreqWin % 1->FG 2->SG 3->LF
    pacValsAllWin = pacValsPerDistPerWin{1,nWin};
    
    %************* Reduction Percentage *****************
    pacValsRaw = squeeze(pacValsAllWin(1,:,:));
    pacValsAfterMP = squeeze(pacValsAllWin(2,:,:));
    pacVals400 = squeeze(pacValsAllWin(3,:,:));

    % before and after MP
    pacValsWinDiffMP = pacValsRaw -  pacValsAfterMP;
    reductionPercentMP = 100*mean(pacValsWinDiffMP./pacValsRaw,'all');
    
    % before and after 400 Elec distance
    pacValsWinDiff400 = pacValsRaw -  pacVals400;
    reductionPercent400 = 100*mean(pacValsWinDiff400./pacValsRaw,'all');

    %************ Statistical analysis *****************
    % before and after MP
    [~,pvalMP,~,statsMP] = ttest(pacValsRaw(:),pacValsAfterMP(:));
    [hMP, ~, ~, ~]=fdr_bh(pvalMP,0.05,'pdep','yes');

    % before and after 400 Elec distance
    [~,pval400,~,stats400] = ttest2(pacValsRaw(:),pacVals400(:));
    [h400, ~, ~, ~]=fdr_bh(pval400,0.05,'pdep','yes');

    %************* metrics table *******************
    metricsMPM1(nWin,:) = [reductionPercentMP statsMP.tstat pvalMP hMP];
    metrics400M1(nWin,:) = [reductionPercent400 stats400.tstat pval400 h400];

end

% Violin Plots
clear g
nWin = 1; % Fast gamma
g(3,3)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin}, 'color',pacLabelsPerWin{1,nWin});
g(3,3).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,3).stat_boxplot('width',0.15);
g(3,3).set_color_options('map',colorVals);
g(3,3).set_names('x',[],'y','KLMI','color','Origin');
g(3,3).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,3).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,3).axe_property('YLim',[-0.5e-4 15e-4],'YTick',[0 3e-4 6e-4 9e-4 12e-4]);
g(3,3).set_layout_options('Position',axisPos{nWin},'legend',false);

nWin = 2; % Slow Gamma
g(3,2)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,2).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,2).stat_boxplot('width',0.15);
g(3,2).set_color_options('map',colorVals);
g(3,2).set_names('x',[],'y','KLMI','color','Origin');
g(3,2).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,2).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,2).axe_property('YLim',[-0.5e-4 15e-4],'YTick',[0 3e-4 6e-4 9e-4 12e-4]);
g(3,2).set_layout_options('Position',axisPos{nWin},'legend',false);

nWin = 3; % Low Frequency
g(3,1)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,1).stat_boxplot('width',0.15);
g(3,1).set_color_options('map',colorVals);
g(3,1).set_names('x',[],'y','KLMI','color','Origin');
g(3,1).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,1).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,1).axe_property('YLim',[-0.5e-4 15e-4],'YTick',[0 3e-4 6e-4 9e-4 12e-4]);
g(3,1).set_layout_options('Position',axisPos{nWin},'legend',false);

%% Monkey 2

monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; % SFOri protocol

% Plot Handles
pacPlotsM2 = getPlotHandles(2,3,[0.55 0.4 0.40 0.53],0.005,0);
pacPlotHand = pacPlotsM2;

% load PAC
for i = 1:length(electrodeDistanceVal)
    clear tmpData
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal{i} '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag{i}) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
    
    % 150 to 500Hz 
    pcolor(pacPlotHand(1,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(freqIdx:end),squeeze(mean(tmpData.pac(:,2,freqIdx:end,:),1))); 
    shading(pacPlotHand(1,i),'interp'); box(pacPlotHand(1,i),'off');
    hold(pacPlotHand(1,i), "on");
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,1}(1) stimFreqWin{1,1}(2)], [stimFreqWin{1,1}(3) stimFreqWin{1,1}(4)], stimFreqWinColors(1,:), 2, '-', 'B');
    colormap(pacPlotHand(1,i),'jet');
    clim(pacPlotHand(1,i), [0 climVals1]);
    set(pacPlotHand(1,i),'Xscale','log');
    set(pacPlotHand(1,i),'Yscale','log');
    set(pacPlotHand(1,i),'FontSize',12);

    % 2 to 150Hz
    pcolor(pacPlotHand(2,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq(1:freqIdx),squeeze(mean(tmpData.pac(:,2,1:freqIdx,:),1))); 
    shading(pacPlotHand(2,i),'interp'); box(pacPlotHand(2,i),'off');
    hold(pacPlotHand(2,i), "on");
    makeBox(pacPlotHand(2,i), [stimFreqWin{1,2}(1) stimFreqWin{1,2}(2)], [stimFreqWin{1,2}(3) stimFreqWin{1,2}(4)], stimFreqWinColors(2,:), 2, '-', 'B');
    makeBox(pacPlotHand(2,i), [stimFreqWin{1,3}(1) stimFreqWin{1,3}(2)], [stimFreqWin{1,3}(3) stimFreqWin{1,3}(4)], stimFreqWinColors(3,:), 2, '-', 'B');
    colormap(pacPlotHand(2,i),'jet');
    clim(pacPlotHand(2,i), [0 climVals2]);
    set(pacPlotHand(2,i),'Xscale','log');
    set(pacPlotHand(2,i),'Yscale','log');
    set(pacPlotHand(2,i),'FontSize',12);
    set(pacPlotHand(2,i),'xTick',[10 30 70]); set(pacPlotHand(2,i),'xTickLabel',[10 30 70]); 
    title(pacPlotHand(1,i), [titleNames{i} '(N=' num2str(size(tmpData.pac,1)) ')'], 'Color', colorVals(i,:), 'FontSize', 12);

    if i==1
        set(pacPlotHand(1,i),'yTick',[200 350 450]); set(pacPlotHand(1,i),'yTickLabel',[ 200 350 450]); 
        set(pacPlotHand(2,i),'yTick',[ 10 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[10 50 100 150]); 

        xlabel(pacPlotHand(2,i),'Phase Frequency (Hz)','Fontsize',12);
        ylabel(pacPlotHand(2,i),'Amplitude Frequency (Hz)','Fontsize',12);
    else
        set(pacPlotHand(1,i),'yTick',[ 250 350 450]); set(pacPlotHand(1,i),'yTickLabel',[]); 
        set(pacPlotHand(2,i),'yTick',[ 10 50 100 150]); set(pacPlotHand(2,i),'yTickLabel',[]); 
    end
end

cb=colorbar(pacPlotHand(1,3));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.065 0 0 -0.18];
cb.Ticks = [0 climVals1];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',10,'Rotation',270);

cb=colorbar(pacPlotHand(2,3));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.0705 0 0 -0.18];
cb.Ticks = [0 climVals2];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',10,'Rotation',270);

set(pacPlotHand,'TickLength',[0.03 0.03]);
set(pacPlotHand(2,1),'xTick',[ 5 10 30 70]); set(pacPlotHand(2,1),'xTickLabel',[ 5 10 30 70]);  
set(pacPlotHand(2,2),'xTick',[ 5 10 30 70]); set(pacPlotHand(2,2),'xTickLabel',[ 5 10 30 70]); 
set(pacPlotHand(2,3),'xTick',[ 5 10 30 70]); set(pacPlotHand(2,3),'xTickLabel',[ 5 10 30 70]); 
set(pacPlotHand(1,1),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,1),'yTickLabel',[5 10 20 30 40 70 100 300 400]); 
set(pacPlotHand(1,2),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,2),'yTickLabel',[]);
set(pacPlotHand(1,3),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,3),'yTickLabel',[]);

xlabel(pacPlotHand(1,1),'Phase Frequency (Hz)','Fontsize',12);
ylabel(pacPlotHand(1,1),'Amplitude Frequency (Hz)','Fontsize',12);
annotation( 'textbox', 'String', 'Monkey 2', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight','Bold','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.72,0.3,0.9,0.7] )

% total number of frequency windows
numFreqWin = length(stimFreqWin) ;

% initialise
pacVals = cell(1, numFreqWin);
pacLabelsPerWin = cell(1, numFreqWin);
pacValsPerDistPerWin = cell(1, numFreqWin);

% axis position
clear axisPos
axisPos = {[0.80 0.02 0.17 0.30], [0.66 0.02 0.17 0.30], [0.51 0.02 0.18 0.30]};

% PAC Plots for each frequency window
for nWin = 1:numFreqWin
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

    pacElec = [];   pacElecTmp = []; pacLabels = []; clear pacLabelsTmp; clear pacPerElecDist;
    for i = 1:length(electrodeDistanceVal)
        tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal{i} '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag{i}) '.mat']), 'pac', 'centerPhaseFreq', 'centerAmpFreq');
        pacPerElec = squeeze(mean(tmpData.pac(:,2, yFreqPos, xFreqPos),[3,4]));

        % PAC for each box
        pacPerElecDist(i,:,:) = squeeze(mean(tmpData.pac(:,2,yFreqPos,xFreqPos),1));
        pacElecTmp = pacPerElec(:);
        for nCh = 1:length(pacElecTmp)
            pacLabelsTmp{nCh,1} = xLabels{i};
        end
        pacLabels = cat(1,pacLabels,pacLabelsTmp);
        pacElec = cat(1,pacElec,pacElecTmp);       
    end
    pacVals{1,nWin} = pacElec;
    pacLabelsPerWin{1,nWin} = pacLabels;
    pacValsPerDistPerWin{1,nWin} = pacPerElecDist;
    
end

%%%%%% Stats and Reduction Percentage %%%%%%%%%%%%%%%%%
for nWin = 1:numFreqWin % 1->FG 2->SG 3->LF
    pacValsAllWin = pacValsPerDistPerWin{1,nWin};
    
    %************* Reduction Percentage *****************
    pacValsRaw = squeeze(pacValsAllWin(1,:,:));
    pacValsAfterMP = squeeze(pacValsAllWin(2,:,:));
    pacVals400 = squeeze(pacValsAllWin(3,:,:));

    % before and after MP
    pacValsWinDiffMP = pacValsRaw -  pacValsAfterMP;
    reductionPercentMP = 100*mean(pacValsWinDiffMP./pacValsRaw,'all');
    
    % before and after 400 Elec distance
    pacValsWinDiff400 = pacValsRaw -  pacVals400;
    reductionPercent400 = 100*mean(pacValsWinDiff400./pacValsRaw,'all');

    %************ Statistical analysis *****************
    % before and after MP
    [~,pvalMP,~,statsMP] = ttest(pacValsRaw(:),pacValsAfterMP(:));
    [hMP, ~, ~, ~]=fdr_bh(pvalMP,0.05,'pdep','yes');

    % before and after 400 Elec distance
    [~,pval400,~,stats400] = ttest2(pacValsRaw(:),pacVals400(:));
    [h400, ~, ~, ~]=fdr_bh(pval400,0.05,'pdep','yes');

    %************* metrics table *******************
    metricsMPM2(nWin,:) = [reductionPercentMP statsMP.tstat pvalMP hMP];
    metrics400M2(nWin,:) = [reductionPercent400 stats400.tstat pval400 h400];

end

% Violin Plots
nWin = 1; % Fast Gamma
g(3,6)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin}, 'color',pacLabelsPerWin{1,nWin});
g(3,6).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,6).stat_boxplot('width',0.15);
g(3,6).set_color_options('map',colorVals);
g(3,6).set_names('x',[],'y','KLMI','color','Origin');
g(3,6).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,6).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,6).axe_property('YLim',[-0.5e-4 15e-4],'YTick',[0 3e-4 6e-4 9e-4 12e-4]);
g(3,6).set_layout_options('Position',axisPos{nWin},'legend',false);

nWin = 2; % Slow Gamma
g(3,5)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,5).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,5).stat_boxplot('width',0.15);
g(3,5).set_color_options('map',colorVals);
g(3,5).set_names('x',[],'y','KLMI','color','Origin');
g(3,5).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,5).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,5).axe_property('YLim',[-0.5e-4 15e-4],'YTick',[0 3e-4 6e-4 9e-4 12e-4]);
g(3,5).set_layout_options('Position',axisPos{nWin},'legend',false);

nWin = 3; % Low Frequency
g(3,4)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,4).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,4).stat_boxplot('width',0.15);
g(3,4).set_color_options('map',colorVals);
g(3,4).set_names('x',[],'y','KLMI','color','Origin');
g(3,4).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,4).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,4).axe_property('YLim',[-0.5e-4 15e-4],'YTick',[0 3e-4 6e-4 9e-4 12e-4]);
g(3,4).set_layout_options('Position',axisPos{nWin},'legend',false);
g.draw();

