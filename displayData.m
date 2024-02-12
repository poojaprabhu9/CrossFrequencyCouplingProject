
function displayData(monkeyName,expDate,protocolName,removeEvokedResponse,tapers,modality,sVarName,sPos,oPos,pacMethod,filterName,nSurrogates,useMPFlag,useCorrectionFlag)

folderSave = 'savedData';
electrodeDistanceValList = [{'0'} {'400'}];
numDistances = length(electrodeDistanceValList);
colorNames = jet(numDistances);
colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allPlots = getPlotHandles(2,5,[0.05 0.05 0.9 0.9],0.05,0.08);
labelName = {'d=0 (N=31)', 'd=0.4'};

for i=1:numDistances
    electrodeDistanceVal = electrodeDistanceValList{i};
    tmpData = load(fullfile(folderSave,[monkeyName expDate protocolName '_removeMean' num2str(removeEvokedResponse) '_Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_d' electrodeDistanceVal '_' sVarName num2str(sPos) '_o' num2str(oPos) '_' pacMethod '_' filterName '_nSur' num2str(nSurrogates) '_MP' num2str(useMPFlag) '.mat']));
    
    for j=1:2 % 1 - Baseline, 2 - stimulus
        if j==1
            titleName='Baseline';
        else
            titleName='Stimulus';
        end
        plot(allPlots(j,1),tmpData.ffcFreq,squeeze(mean(tmpData.ffc(:,j,:),1)),'color',colorNames(i,:),'LineWidth',1.5);
        ylabel(allPlots(j,1),'FFC','FontWeight','bold');
        title(allPlots(j,1),titleName,'FontSize',14);
        hold(allPlots(j,1),'on');

        plot(allPlots(j,2),tmpData.sfcFreq,squeeze(mean(tmpData.sfc(:,j,:),1)),'color',colorNames(i,:),'LineWidth',1.5);
        ylabel(allPlots(j,2),'SFC','FontWeight','bold');
        title(allPlots(j,2),titleName,'FontSize',14);
        hold(allPlots(j,2),'on');               

        plot(allPlots(j,3),1000*tmpData.xsSTA,squeeze(mean(tmpData.staVals(:,j,:),1)),'color',colorNames(i,:),'LineWidth',1.5);
        ylabel(allPlots(j,3),'STA','FontWeight','bold');
        title(allPlots(j,3),titleName,'FontSize',14);
        hold(allPlots(j,3),'on');
        
        if useCorrectionFlag
            [h, ~, ~, ~]=fdr_bh(tmpData.pval,0.05,'pdep','yes');
            pac=tmpData.pac.*h;
        else 
            pac=tmpData.pac;
        end
        pcolor(allPlots(j,3+i),tmpData.centerPhaseFreq,tmpData.centerAmpFreq,squeeze(mean(pac(:,j,:,:),1))); 
        title(allPlots(j,3+i),[titleName ' for ' labelName{i}],'FontSize',14);
        ylabel(allPlots(j,3+i),'Amplitude Frequency (Hz)','FontWeight','bold');
        clim(allPlots(j,3+i),[min(pac(:)) max(pac(:))]);
        shading(allPlots(j,3+i),'interp');
        cb=colorbar(allPlots(j,3+i));
        xlabel(cb,pacMethod,'FontWeight','bold','FontSize',10,'Rotation',270);
        if j==2
            xlabel(allPlots(j,1),'Frequency (Hz)','FontWeight','bold');
            xlabel(allPlots(j,2),'Frequency (Hz)','FontWeight','bold');
            xlabel(allPlots(j,3),'Time (ms)','FontWeight','bold');
            xlabel(allPlots(j,3+i),'Phase Frequency (Hz)','FontWeight','bold');
        end        
        axis(allPlots(j,1),[0 tmpData.ffcFreq(end) 0 1]);
        axis(allPlots(j,2),[0 tmpData.sfcFreq(end) 0 0.4]);
        axis(allPlots(j,3),[tmpData.xsSTA(1)*10^3 tmpData.xsSTA(end)*10^3 -50 10]);
    end
legend(allPlots(1,2),labelName,'fontsize',12 ,'box','off'); 
end

end