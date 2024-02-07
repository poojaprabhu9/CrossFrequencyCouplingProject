
function displayData(monkeyName,expDate,protocolName,removeEvokedResponse,tapers,modality,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag)

folderSave = 'savedData';
electrodeDistanceValList = [{'0'} {'400'}];
numDistances = length(electrodeDistanceValList);
colorNames = jet(numDistances);
colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allPlots = getPlotHandles(2,5,[0.05 0.05 0.9 0.9],0.05,0.05);

for i=1:numDistances
    electrodeDistanceVal = electrodeDistanceValList{i};
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponse) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sPos) '_o' num2str(oPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']));

    for j=1:2 % 1 - Baseline, 2 - stimulus
        plot(allPlots(j,1),tmpData.ffcFreq,squeeze(mean(tmpData.ffc(:,j,:),1)),'color',colorNames(i,:));
        hold(allPlots(j,1),'on');

        plot(allPlots(j,2),tmpData.sfcFreq,squeeze(mean(tmpData.sfc(:,j,:),1)),'color',colorNames(i,:));
        hold(allPlots(j,2),'on');

        plot(allPlots(j,3),1000*tmpData.xsSTA,squeeze(mean(tmpData.staVals(:,j,:),1)),'color',colorNames(i,:));
        hold(allPlots(j,3),'on');

        pcolor(allPlots(j,4+i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq,squeeze(mean(tmpData.pac(:,j,:,:),1))); 
        shading(allPlots(j,4+i),'interp');

    end
end