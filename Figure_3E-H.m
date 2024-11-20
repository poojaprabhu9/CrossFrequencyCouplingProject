% Figure 3E-H: Plotting for Phase-Amplitude Coupling (PAC) for 'all-all'
% conditions for Electrocorticography (ECoG), under two cases,
% 1. D0 (single electrode LFP)
% 2. D0 after MP (single electrode LFP after Matching Pursuit(MP))
% For the three cases the violin plots are generated for three different
% frequency bands (Low Frequency, Slow Gamma, and Fast Gamma)

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms
% gramm data visualization toolbox: https://github.com/piermorel/gramm

clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Electrode choices
modality = 'ECoG'; 

% Datapath
folderSave = 'savedData';

% Colormap for PAC Plots
colormap jet

% SF and Ori position for 'all-all' conditions
sfPos = 6; oriPos = 9;

% Gamma range
gammaRange = {[40 70], [25 40]}; % Fast Gamma and Slow Gamma

% Frequency range for three bounding boxes-Fast gamma, slow gamma and low
% frequency, respectively
stimFreqWin = {[gammaRange{2}(1) gammaRange{2}(2) 87 157], [4 8 57 177] };

% Color for two boxes- slow gamma and low frequency, respectively
stimFreqWinColors = [92/255 64/255 51/255; 0 0 0];

% labels for each frequency windows
labelNames = { {'Slow Gamma','with 80-150 Hz'}, {'Low Frequency(<10 Hz)','with 50-180 Hz'}};

% labels needed to create the structure
xLabels = {'D0', 'D0_afterMP'};

% color vals for two cases (D0, D0 after MP)
colorVals = [1 0 1; 0.2 0.8 0.5]; 

% Signal processing choices
removeEvokedResponseFlag = 1; 
tapers = [1 1];
pacMethod = 'klmi'; 
filterName = 'fir'; 
nSurrogates = 0; % Number of surrogates
useMPFlag = {'0','1'}; % Particular to MP filtered data
sVarName = 'sf';
electrodeDistanceVal = {'0','0'};
titleNames = {'D_{0}', 'D_{0} after MP'};

%% Monkey 1

monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol

% Plot Handles
pacPlotsM1 = getPlotHandles(1,2,[0.05 0.4 0.40 0.5],0.005,0);
pacPlotHand = pacPlotsM1;

% Color bar limits for 2-500 Hz 
climVals1 = 0.001; 

% load PAC
for i = 1:length(electrodeDistanceVal)
    clear tmpData
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal{i} '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag{i}) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
    pcolor(pacPlotHand(1,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq,squeeze(mean(tmpData.pac(:,2,:,:),1))); 
    shading(pacPlotHand(1,i),'interp'); box(pacPlotHand(1,i),'off');
    hold(pacPlotHand(1,i), "on");
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,1}(1) stimFreqWin{1,1}(2)], [stimFreqWin{1,1}(3) stimFreqWin{1,1}(4)], stimFreqWinColors(1,:), 2, '-', 'B');                        
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,2}(1) stimFreqWin{1,2}(2)], [stimFreqWin{1,2}(3) stimFreqWin{1,2}(4)], stimFreqWinColors(2,:), 2, '-', 'B');
    colormap(pacPlotHand(1,i),'jet');
    clim(pacPlotHand(1,i), [0 climVals1]);
    set(pacPlotHand(1,i),'Xscale','log');
    set(pacPlotHand(1,i),'Yscale','log');
    title(pacPlotHand(1,i), [titleNames{i} '(N=' num2str(size(tmpData.pac,1)) ')'], 'Color', colorVals(i,:), 'FontSize', 20);
end

cb=colorbar(pacPlotHand(1,2));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.07 0 0 -0.4];
cb.Ticks = [0 climVals1];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',10,'Rotation',270);

set(pacPlotHand,'TickLength',[0.03 0.03]);
set(pacPlotHand(1,1),'xTick',[ 5 10 30 70]); set(pacPlotHand(1,1),'xTickLabel',[ 5 10 30 70]);  
set(pacPlotHand(1,2),'xTick',[ 5 10 30 70]); set(pacPlotHand(1,2),'xTickLabel',[ 5 10 30 70]); 

set(pacPlotHand(1,1),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,1),'yTickLabel',[5 10 20 30 40 70 100 300 400]); 
set(pacPlotHand(1,2),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,2),'yTickLabel',[]);

set(pacPlotHand(1,1),'FontSize',12); set(pacPlotHand(1,2),'FontSize',12);  
xlabel(pacPlotHand(1,1),'Phase Frequency (Hz)','Fontsize',12);
ylabel(pacPlotHand(1,1),'Amplitude Frequency (Hz)','Fontsize',12);
annotation( 'textbox', 'String', 'Monkey 1', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight','Bold', 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.22,0.3,0.9,0.7] )

% total number of frequency windows
numFreqWin = length(stimFreqWin) ;

% initialise
pacVals = cell(1, numFreqWin);
pacLabelsPerWin = cell(1, numFreqWin);
pacValsPerDistPerWin = cell(1, numFreqWin);

% axis position
axisPos = {[0.26 0.02 0.17 0.30], [0.06 0.02 0.17 0.30]};

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
    pacElec = [];   pacElecTmp = []; pacLabels = []; clear pacLabelsTmp; clear pacPerElecDist;
    for i = 1:length(electrodeDistanceVal)
        tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal{i} '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag{i}) '.mat']), 'pac', 'centerPhaseFreq', 'centerAmpFreq');
        xFreqPos = intersect(find(tmpData.centerPhaseFreq>= xFreq1),find(tmpData.centerPhaseFreq<=  xFreq2));
        yFreqPos = intersect(find(tmpData.centerAmpFreq>= stimFreqWin{1,nWin}(3)),find(tmpData.centerAmpFreq<= stimFreqWin{1,nWin}(4)));

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
for nWin = 1:numFreqWin % 1->SG 2->LF
    pacValsAllWin = pacValsPerDistPerWin{1,nWin};
    
    %************* Reduction Percentage *****************
    pacValsRaw = squeeze(pacValsAllWin(1,:,:));
    pacValsAfterMP = squeeze(pacValsAllWin(2,:,:));

    % before and after MP
    pacValsWinDiffMP = pacValsRaw -  pacValsAfterMP;
    reductionPercentMP = 100*mean(pacValsWinDiffMP./pacValsRaw,'all');
    
    %************ Statistical analysis *****************
    % before and after MP
    [~,pvalMP,~,statsMP] = ttest(pacValsRaw(:),pacValsAfterMP(:));
    [hMP, ~, ~, ~]=fdr_bh(pvalMP,0.05,'pdep','yes');

    %************* metrics table *******************
    metricsMPM1(nWin,:) = [reductionPercentMP statsMP.tstat pvalMP hMP];

end

% Violin Plots
clear g
nWin = 1; % Slow Gamma
g(3,2)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,2).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,2).stat_boxplot('width',0.15);
g(3,2).set_color_options('map',colorVals);
g(3,2).set_names('x',[],'y','KLMI','color','Origin');
g(3,2).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,2).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,2).axe_property('YLim',[-0.5e-4 25e-4],'YTick',[0 5e-4 10e-4 15e-4 20e-4]);
g(3,2).set_layout_options('Position',axisPos{nWin},'legend',false);

nWin = 2; % Low frequency
g(3,1)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,1).stat_boxplot('width',0.15);
g(3,1).set_color_options('map',colorVals);
g(3,1).set_names('x',[],'y','KLMI','color','Origin');
g(3,1).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,1).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,1).axe_property('YLim',[-0.5e-4 25e-4],'YTick',[0 5e-4 10e-4 15e-4 20e-4]);
g(3,1).set_layout_options('Position',axisPos{nWin},'legend',false);

%% Monkey 2

monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; % SF-Ori protocol

% Plot Handles
pacPlotsM2 = getPlotHandles(1,2,[0.55 0.4 0.40 0.5],0.005,0);
pacPlotHand = pacPlotsM2;

% Colorbar limits for PAC plots
climVals1 = 0.001; 

% load PAC
for i = 1:length(electrodeDistanceVal)
    clear tmpData
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponseFlag) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal{i} '_' sVarName num2str(sfPos) '_o' num2str(oriPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag{i}) '.mat']), 'pac', 'centerAmpFreq', 'centerPhaseFreq');
    pcolor(pacPlotHand(1,i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq,squeeze(mean(tmpData.pac(:,2,:,:),1))); 
    shading(pacPlotHand(1,i),'interp'); box(pacPlotHand(1,i),'off');
    hold(pacPlotHand(1,i), "on");
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,1}(1) stimFreqWin{1,1}(2)], [stimFreqWin{1,1}(3) stimFreqWin{1,1}(4)], stimFreqWinColors(1,:), 2, '-', 'B');                        
    makeBox(pacPlotHand(1,i), [stimFreqWin{1,2}(1) stimFreqWin{1,2}(2)], [stimFreqWin{1,2}(3) stimFreqWin{1,2}(4)], stimFreqWinColors(2,:), 2, '-', 'B');
    colormap(pacPlotHand(1,i),'jet');
    clim(pacPlotHand(1,i), [0 climVals1]);
    set(pacPlotHand(1,i),'Xscale','log');
    set(pacPlotHand(1,i),'Yscale','log');
    title(pacPlotHand(1,i), [titleNames{i} '(N=' num2str(size(tmpData.pac,1)) ')'], 'Color', colorVals(i,:), 'FontSize', 20);
end

cb=colorbar(pacPlotHand(1,2));
cb.Ruler.Exponent = -3;
cb.Position = cb.Position + [0.07 0 0 -0.4];
cb.Ticks = [0 climVals1];
cb.TickDirection = 'none';
xlabel(cb,'KLMI','FontWeight','bold','FontSize',10,'Rotation',270);

set(pacPlotHand,'TickLength',[0.03 0.03]);
set(pacPlotHand(1,1),'xTick',[ 5 10 30 70]); set(pacPlotHand(1,1),'xTickLabel',[ 5 10 30 70]);  
set(pacPlotHand(1,2),'xTick',[ 5 10 30 70]); set(pacPlotHand(1,2),'xTickLabel',[ 5 10 30 70]); 

set(pacPlotHand(1,1),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,1),'yTickLabel',[5 10 20 30 40 70 100 300 400]); 
set(pacPlotHand(1,2),'yTick',[5 10 20 30 40 70 100 300 400]); set(pacPlotHand(1,2),'yTickLabel',[]);
 
set(pacPlotHand(1,1),'FontSize',12); set(pacPlotHand(1,2),'FontSize',12);  
xlabel(pacPlotHand(1,1),'Phase Frequency (Hz)','Fontsize',12);
ylabel(pacPlotHand(1,1),'Amplitude Frequency (Hz)','Fontsize',12);
annotation( 'textbox', 'String', 'Monkey 2', 'Color', 'black', ...
            'FontSize', 14, 'FontWeight','Bold','Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.73,0.3,0.9,0.7] )

% total number of frequency windows
numFreqWin = length(stimFreqWin) ;

% initialise
pacVals = cell(1, numFreqWin);
pacLabelsPerWin = cell(1, numFreqWin);
pacValsPerDistPerWin = cell(1, numFreqWin);

clear axisPos
axisPos = {[0.76 0.02 0.17 0.30], [0.55 0.02 0.17 0.30]};

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
for nWin = 1:numFreqWin % 1->SG 2->LF
    pacValsAllWin = pacValsPerDistPerWin{1,nWin};
    
    %************* Reduction Percentage *****************
    pacValsRaw = squeeze(pacValsAllWin(1,:,:));
    pacValsAfterMP = squeeze(pacValsAllWin(2,:,:));

    % before and after MP
    pacValsWinDiffMP = pacValsRaw -  pacValsAfterMP;
    reductionPercentMP = 100*mean(pacValsWinDiffMP./pacValsRaw,'all');
    
    %************ Statistical analysis *****************
    % before and after MP
    [~,pvalMP,~,statsMP] = ttest(pacValsRaw(:),pacValsAfterMP(:));
    [hMP, ~, ~, ~]=fdr_bh(pvalMP,0.05,'pdep','yes');

    %************* metrics table *******************
    metricsMPM2(nWin,:) = [reductionPercentMP statsMP.tstat pvalMP hMP];

end

% Violin Plots
nWin = 1; % Slow Gamma
g(3,5)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,5).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,5).stat_boxplot('width',0.15);
g(3,5).set_color_options('map',colorVals);
g(3,5).set_names('x',[],'y','KLMI','color','Origin');
g(3,5).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,5).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,5).axe_property('YLim',[-0.5e-4 20e-4],'YTick',[0 5e-4 10e-4 15e-4]);
g(3,5).set_layout_options('Position',axisPos{nWin},'legend',false);

nWin = 2; % Low frequency
g(3,4)=gramm('x',pacLabelsPerWin{1,nWin},'y',pacVals{1,nWin},'color',pacLabelsPerWin{1,nWin});
g(3,4).stat_violin('normalization','area','dodge',0,'fill','edge');
g(3,4).stat_boxplot('width',0.15);
g(3,4).set_color_options('map',colorVals);
g(3,4).set_names('x',[],'y','KLMI','color','Origin');
g(3,4).set_title(labelNames{nWin},'color',stimFreqWinColors(nWin,:),'FontSize',12);
g(3,4).axe_property('XTick',[],'XTickLabel',[],'Ygrid','on');
g(3,4).axe_property('YLim',[-0.5e-4 20e-4],'YTick',[0 5e-4 10e-4 15e-4]);
g(3,4).set_layout_options('Position',axisPos{nWin},'legend',false);
g.draw();




