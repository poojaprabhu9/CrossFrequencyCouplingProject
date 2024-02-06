
function displayData(monkeyName,expDate,protocolName,removeEvokedResponse,tapers,modality,sVarName,sPos,oPos,methodVar,filterMethod,nSurr,useMPFlag,metrics,correction,pThresh)

folderSave = 'savedData';
electrodeDistanceValList = [{'0'} {'400'}];
numDistances = length(electrodeDistanceValList);
colorNames = jet(numDistances);
colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allPlots = getPlotHandles(2,4,[0.05 0.05 0.9 0.9],0.05,0.05);

for i=1:numDistances
    electrodeDistanceVal = electrodeDistanceValList{i};
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponse) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sPos) '_o' num2str(oPos) '_' methodVar '_' filterMethod '_nSur' num2str(nSurr) '_MP' num2str(useMPFlag) '.mat']));

    for j=1:2 % 1 - Baseline, 2 - stimulus
        plot(allPlots(j,1),tmpData.ffcFreq,squeeze(mean(tmpData.ffc(:,j,:),1)),'color',colorNames(i,:));
        hold(allPlots(j,1),'on');

        plot(allPlots(j,2),1000*tmpData.xsSTA,squeeze(mean(tmpData.staVals(:,j,:),1)),'color',colorNames(i,:));
        hold(allPlots(j,2),'on');

        pcolor(allPlots(j,2+i),tmpData.freqVecPhase,tmpData.freqVecAmp,squeeze(mean(tmpData.pacmat(:,j,:,:),1))); 
        shading(allPlots(j,2+i),'interp');
        colorbar(allPlots(j,2+i));
    end
end