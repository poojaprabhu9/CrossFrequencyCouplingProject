
function displayData(monkeyName,expDate,protocolName,removeEvokedResponse,tapers,spikeElectrodes,analysisPeriodNum,fPos)

if ~exist('removeEvokedResponse','var'); removeEvokedResponse=1;        end
if ~exist('tapers','var'); tapers=[1 1];                                end
if ~exist('spikeElectrodes','var'); spikeElectrodes=[];                 end
if ~exist('analysisPeriodNum','var'); analysisPeriodNum=2;              end
if ~exist('fPos','var'); fPos=6;                                        end

rfData = load([monkeyName 'MicroelectrodeRFData.mat']);

goodElectrodes = rfData.highRMSElectrodes;
lfpElectrodes = goodElectrodes(goodElectrodes<=81);

if isempty(spikeElectrodes)
    spikeElectrodes = lfpElectrodes;
end

folderSave = 'savedData';

electrodeDistanceList{1} = 0;
electrodeDistanceList{2} = 400;
electrodeDistanceList{3} = [400 1600];
electrodeDistanceList{4} = [1600 2400];
electrodeDistanceList{5} = [2400 4000];

distanceColorNames = ['k';'r';'b';'c';'g'];
numCombinations = length(electrodeDistanceList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
numElectrodes = length(spikeElectrodes);

numOrientations = 9;
plotHandlesFFC = getPlotHandles(1,numOrientations,[0.05 0.85 0.9 0.125],0.01); linkaxes(plotHandlesFFC);
plotHandlesSFC = getPlotHandles(1,numOrientations,[0.05 0.7 0.9 0.125],0.01); linkaxes(plotHandlesSFC);
plotHandlesSTA = getPlotHandles(1,numOrientations,[0.05 0.55 0.9 0.1],0.01); linkaxes(plotHandlesSTA);

ffcData = cell(1,numCombinations);
sfcData = cell(1,numCombinations);
staData = cell(1,numCombinations);

for i=1:numElectrodes
    spikeElectrode = spikeElectrodes(i);

    disp(['Working on electrode: ' num2str(spikeElectrode)]);
    tmp = load(fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(spikeElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '.mat']));

    electrodesToCombine = getElectrodeDistance(lfpElectrodes,spikeElectrode,electrodeDistanceList);

    for j=1:numCombinations
        ffcData{j} = cat(1,ffcData{j},(tmp.ffc(electrodesToCombine{j},fPos,:,analysisPeriodNum,:)));
        sfcData{j} = cat(1,sfcData{j},(tmp.sfc(electrodesToCombine{j},fPos,:,analysisPeriodNum,:)));
        staData{j} = cat(1,staData{j},cell2mat((tmp.staVals(electrodesToCombine{j},fPos,:,analysisPeriodNum))));
    end
end

ffcFreq = tmp.ffcFreq;
sfcFreq = tmp.sfcFreq;
xsSTA = 1000*tmp.xsSTA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for o=1:numOrientations
    for i=1:numCombinations
        
        plot(plotHandlesFFC(o),ffcFreq,squeeze(mean(ffcData{i}(:,:,o,:,:),1)),'color',distanceColorNames(i));
        hold(plotHandlesFFC(o),'on');

        plot(plotHandlesSFC(o),sfcFreq,squeeze(mean(sfcData{i}(:,:,o,:,:),1)),'color',distanceColorNames(i));
        hold(plotHandlesSFC(o),'on');

        plot(plotHandlesSTA(o),xsSTA,squeeze(mean(staData{i}(:,:,o,:,:),1)),'color',distanceColorNames(i));
        hold(plotHandlesSTA(o),'on');
    end

    set(plotHandlesFFC(o),'XTickLabel',[]);
    if o==1
        ylabel(plotHandlesFFC(o),'FFC');
        ylabel(plotHandlesSFC(o),'SFC');
        ylabel(plotHandlesSTA(o),'STA');
    else
        set(plotHandlesFFC(o),'YTickLabel',[]);
        set(plotHandlesSFC(o),'YTickLabel',[]);
        set(plotHandlesSTA(o),'YTickLabel',[]);
    end

    xlabel(plotHandlesSFC(o),'Frequency (Hz)');
    xlabel(plotHandlesSTA(o),'Time (ms)');

    axis(plotHandlesFFC(o),[0 ffcFreq(end) 0 1]);
    axis(plotHandlesSFC(o),[0 ffcFreq(end) 0 0.2]);
    axis(plotHandlesSTA(o),[xsSTA(1) xsSTA(end) -40 10]);
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