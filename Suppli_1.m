% Suppli 1A or B: Time Frequency (TF) Plots for all the spatial
% frequencies and orientations
% Script plots the TF plots for one monkey at a time 

% Need to have the following folders in Matlab's path
% CommonPrograms: https://github.com/supratimray/CommonPrograms

clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%% Choice of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index values of all the spatial frequencies and Orientations
sfPos = [1 2 3 4 5 6];
oriPos = [1 2 3 4 5 6 7 8 9];

% Label names
sfVals = {'0.5', '1', '2', '4', '8', 'all'};
oriVals = {'0', '22.5', '45', '67.5', '90', '112.5', '135', '157.5', 'all'};

% Electrode choices
modality = 'LFP'; %'LFP' or 'ECoG'

% Signal processing choices
tapers = [2 3];

% Frequency ranges
gammaRange = {[40 70], [25 40]}; % Fast Gamma and Slow Gamma; Murty et al., 2018

% Datapath
folderSave = 'savedData';

% Time limits 
timeLim = [-0.2220 1.0030];

% Frequency Limits
freqLim = [0 100];

% Baseline period
blRange=[-0.5 0];

% colormap for PAC Plots
colormap jet

% Protocol details
% comment the below line according to chosen monkey
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002'; % SFOri protocol
monkeyName = 'kesariH'; expDate = '270218'; protocolName = 'GRF_001'; 

% Plot Handles
figure(1)
tfPlots = getPlotHandles(6,9,[0.09 0.08 0.85 0.85],0.01,0.03);

% load PSD data
if strcmp(modality,'LFP')
    tfData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_spike' modality '_TimeFreq_For_All_Conditions.mat']));
elseif strcmp(modality,'ECoG')
    tfData = load(fullfile(folderSave,[monkeyName expDate protocolName '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_TimeFreq_For_All_Conditions.mat']));
end
tf = tfData.energyValues;
tfFreqVals = tfData.mtFreq;
tfTimeVals = tfData.mtTime ;

freqPos = intersect(find(tfFreqVals>= freqLim(1)),find(tfFreqVals<=freqLim(2)));
timePos = intersect(find(tfTimeVals >= timeLim(1)),find(tfTimeVals<= timeLim(2)));

for sPos = 1:length(sfPos)
    for oPos = 1:length(oriPos)
        tfValsPerElec = zeros(size(tf,1),length(tfTimeVals),length(tfFreqVals));
        for nCh = 1:size(tf,1)
            tmpPower=squeeze(tf(nCh,sfPos(sPos),oriPos(oPos),:,:));

            blPos = intersect(find(tfTimeVals>=blRange(1)),find(tfTimeVals<blRange(2)));
            logTmpPower = log10(tmpPower);
            blPower = mean(logTmpPower(blPos,:),1);
            logSBL = repmat(blPower,length(tfTimeVals),1); %#ok<NASGU>
            deltaTF = 10*(logTmpPower-logSBL);
            tfValsPerElec(nCh,:,:)=deltaTF;
        end
        avgTfValsPerElec=squeeze(mean(tfValsPerElec,1));
        pcolor(tfPlots(sPos,oPos),tfTimeVals(timePos),tfFreqVals(freqPos),avgTfValsPerElec(timePos,freqPos)'); 

       shading(tfPlots(sPos,oPos),'interp'); box(tfPlots(sPos,oPos),'off');
       colormap(tfPlots(sPos,oPos),'jet');

       % comment the below line according to chosen monkey
       clim(tfPlots(sPos,oPos), [-3 7]); % for Monkey 1 (alpaH)
       clim(tfPlots(sPos,oPos), [-2 7]); % for Monkey 2 (kesariH)

       % Plot stim time range
       xline(tfPlots(sPos,oPos),0.2530,'Color',[0 0 0],'LineWidth',0.5)
       xline(tfPlots(sPos,oPos),0.7280,'Color',[0 0 0],'LineWidth',0.5)

       % Plot gamma Range
       % fast gamma
       yline(tfPlots(sPos,oPos),gammaRange{1}(1),'Color',[153/255 102/255 204/255],'LineWidth',1);yline(tfPlots(sPos,oPos),gammaRange{1}(2),'Color',[153/255 102/255 204/255],'LineWidth',1);
       % slow gamma
       yline(tfPlots(sPos,oPos),gammaRange{2}(1),'Color',[92/255 64/255 51/255],'LineWidth',1);yline(tfPlots(sPos,oPos),gammaRange{2}(2),'Color',[92/255 64/255 51/255],'LineWidth',1);

        if sPos == 6
            xlabel(tfPlots(sPos,oPos), 'Time (s)','FontSize',12)
            set(tfPlots(sPos,oPos),'FontSize',12);
        else
            set(tfPlots(sPos,oPos),'xTickLabel',[]);   
        end
        if oPos == 1
            set(tfPlots(sPos,oPos),'yTick',[0 50 100]); 
            set(tfPlots(sPos,oPos),'yTickLabel',[0 50 100]); 
            set(tfPlots(sPos,oPos),'FontSize',12);
        else
            set(tfPlots(sPos,oPos),'yTick',[0 50 100]); 
            set(tfPlots(sPos,oPos),'yTickLabel',[]); 
        end
        if sPos==1
           title(tfPlots(sPos,oPos), oriVals{oPos}, 'FontSize', 12);
        end
    end
end
ylabel(tfPlots(6,1), 'Frequency (Hz)','FontSize', 12)
cb=colorbar(tfPlots(6,9));
cb.Position = cb.Position + [0.065 0 0 -0.07];
% comment the below line according to chosen monkey
cb.Ticks = [-3 7]; % for Monkey 1 (alpaH)
cb.Ticks = [-2 7]; % for Monkey 2 (kesariH)
cb.TickDirection = 'none';
for j = 1:9
    set(tfPlots(6,j),'xTickLabel',[]);
    set(tfPlots(6,j),'xTick',[0 0.4 0.8]);
    set(tfPlots(6,j),'xTickLabel',[0 0.4 0.8]);
end
annotation( 'textbox', 'String', 'Change in power(dB)', 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.942,0.14,0.2,0], 'rotation',90, 'FontWeight','bold')

annotation( 'textbox', 'String', sfVals{1}, 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.035,0.9,0.9,0], 'rotation',0, 'FontWeight','bold')
annotation( 'textbox', 'String', sfVals{2}, 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.035,0.78,0.9,0], 'rotation',0, 'FontWeight','bold')
annotation( 'textbox', 'String', sfVals{3}, 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.035,0.62,0.9,0], 'rotation',0, 'FontWeight','bold')
annotation( 'textbox', 'String', sfVals{4}, 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.035,0.48,0.9,0], 'rotation',0, 'FontWeight','bold')
annotation( 'textbox', 'String', sfVals{5}, 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.035,0.32,0.9,0], 'rotation',0, 'FontWeight','bold')
annotation( 'textbox', 'String', sfVals{6}, 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.035,0.20,0.9,0], 'rotation',0, 'FontWeight','bold')

annotation( 'textbox', 'String', 'Spatial Frequency (cpd)', 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.02,0.40,0.9,0], 'rotation',90, 'FontWeight','bold')

annotation( 'textbox', 'String', 'Orientation (degrees)', 'Color', 'black', ...
            'FontSize', 12, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.45,0.10,0.9,0.9], 'rotation',0, 'FontWeight','bold')


