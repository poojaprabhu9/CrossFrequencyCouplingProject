
function displayData(monkeyName,expDate,protocolName,removeEvokedResponse,tapers,spikeElectrodes,analysisPeriodNum,modality, electrodeDistanceVal, tuning, filterMethod, methodVar, metrics, correction, pThresh)

if ~exist('removeEvokedResponse','var'); removeEvokedResponse=1;        end
if ~exist('tapers','var'); tapers=[1 1];                                end
if ~exist('spikeElectrodes','var'); spikeElectrodes=[];                 end
if ~exist('analysisPeriodNum','var'); analysisPeriodNum=2;              end
if ~exist('tuning','var'); tuning='Ori';                                end
if ~exist('electrodeDistanceVal','var'); electrodeDistanceVal='0';      end


if (strcmp(tuning,'SF4Ori90') || strcmp(tuning,'SF2Ori90'))
    condition='single';
    if strcmp(metrics,'others')
        electrodeDistanceVal={'0','400'};
    end
    numCombinations=length(electrodeDistanceVal);
else
    condition='multiple';
    numCombinations=1;
end

% load highRMS electrodes
rfData = load([monkeyName 'MicroelectrodeRFData.mat']);

goodElectrodes = rfData.highRMSElectrodes;

% Choice of good Electrodes based on chosen modality
switch modality

    case 'ECoG'

        if strcmp(monkeyName,'alpaH')

            ecogElectrodes=[82 85 86 88 89];

        elseif strcmp(monkeyName,'kesariH')

            ecogElectrodes=[85 86 88 89];
        end
        lfpElectrodes=ecogElectrodes;
        spikeElectrodes=ecogElectrodes;

    case 'LFP'

        lfpElectrodes = goodElectrodes(goodElectrodes<=81);

end

if isempty(spikeElectrodes)
    spikeElectrodes = lfpElectrodes;
end

folderSave = 'savedData';

distanceColorNames = ['k';'r';'b';'c';'g'];

numCombinations=length(electrodeDistanceVal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
% Based on type of conditions, Figures are generated
switch condition
    case 'single'  % SF4Ori90 or SF2Ori90      

        if strcmp(metrics,'others') % plotting FFC, SFC, STA
            labelNames={'LFP-bsl','LFP-stim','ECoG-bsl','ECoG-stim'};
            numOrientations=4; % LFP-bsl,LFP-stim,ECoG-bsl,ECoG-stim,

            modality={'LFP','ECoG'};
            analysisPeriodNum=[1 2];

            plotHandlesFFC = getPlotHandles(1,numOrientations,[0.05 0.85 0.9 0.125],0.01); linkaxes(plotHandlesFFC);
            plotHandlesSFC = getPlotHandles(1,numOrientations,[0.05 0.70 0.9 0.125],0.01); linkaxes(plotHandlesSFC);
            plotHandlesSTA = getPlotHandles(1,numOrientations,[0.05 0.55 0.9 0.1],0.01); linkaxes(plotHandlesSTA);
            ffcData = cell(1,numCombinations);
            sfcData = cell(1,numCombinations);
            staData = cell(1,numCombinations);
            

            for moda=1:length(modality)

                switch modality{moda}
    
                        case 'ECoG'
                    
                            if strcmp(monkeyName,'alpaH')
                    
                                ecogElectrodes=[82 85 86 88 89];
                    
                            elseif strcmp(monkeyName,'kesariH')
                    
                                ecogElectrodes=[85 86 88 89];
                            end
                            lfpElectrodes=ecogElectrodes;
                            
                            electrodeDistanceVal={'0','400'};
                            numCombinations=length(electrodeDistanceVal);
    
                            numElectrodes=length(ecogElectrodes);
                            for dist=1:numCombinations
                                if strcmp(electrodeDistanceVal{dist},'400')
                                    ffcData{dist}(3,:)=0; %average across all ecog electrodes in baseline
                                    ffcData{dist}(4,:)=0; %average across all ecog electrodes in stim
                
                                    sfcData{dist}(3,:)=0; %average across all ecog electrodes in baseline
                                    sfcData{dist}(4,:)=0; %average across all ecog electrodes in stim
                
                                    staData{dist}(3,:)=0; %average across all ecog electrodes in baseline
                                    staData{dist}(4,:)=0; %average across all ecog electrodes in stim



                                else
                                    for n=1:numElectrodes
                                        ecogElectrode=ecogElectrodes(n);
            
                                        disp(['Working on electrode: ' num2str(ecogElectrode)]);
                                        tmp = load(fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(ecogElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality{moda} '_' tuning 'tuning_Allmethods_for_distance' electrodeDistanceVal{dist} '.mat']));
                                        blah=cell2mat(tmp.staVals);
                                    
                                      
                                    
                                       for ee=1:size(blah,1)
                                            for tt=1:size(blah,4)                       
                                                staValsMatEcog(ee,1,1,tt,:)=squeeze(blah(ee,:,1,tt));
                                            end
                                       end
            
                                       for ee=1:size(tmp.ffc,1) % number of LFP elec
                                        for tt=1:size(tmp.ffc,4) % number of analysis period                           
                            
                                            ecog_ffc(n,ee,1,1,tt,:)=tmp.ffc(ee,1,1,tt,:);
                                            ecog_sfc(n,ee,1,1,tt,:)=tmp.sfc(ee,1,1,tt,:);
                                            ecog_staValsMat(n,ee,1,1,tt,:)=staValsMatEcog(ee,1,1,tt,:);
                                        end
                                       end
                                    end
    
                                    ffcData{dist}(3,:)=squeeze(mean(ecog_ffc(:,:,1,1,1,:),[1 2])); %average across all ecog electrodes in baseline
                                    ffcData{dist}(4,:)=squeeze(mean(ecog_ffc(:,:,1,1,2,:),[1 2])); %average across all ecog electrodes in stim
                    
                                    sfcData{dist}(3,:)=squeeze(mean(ecog_sfc(:,:,1,1,1,:),[1 2])); %average across all ecog electrodes in baseline
                                    sfcData{dist}(4,:)=squeeze(mean(ecog_sfc(:,:,1,1,2,:),[1 2])); %average across all ecog electrodes in stim
                    
                                    staData{dist}(3,:)=squeeze(mean(ecog_staValsMat(:,:,1,1,1,:),[1 2])); %average across all ecog electrodes in baseline
                                    staData{dist}(4,:)=squeeze(mean(ecog_staValsMat(:,:,1,1,2,:),[1 2])); %average across all ecog electrodes in stim
    
                                    end
                            end
                    
                        case 'LFP'
                    
                            lfpElectrodes = goodElectrodes(goodElectrodes<=81);
                            spikeElectrodes=spikeElectrodes;
                            electrodeDistanceVal={'0','400'};
                            numCombinations=length(electrodeDistanceVal);
                            numElectrodes=length(spikeElectrodes);
                            for dist=1:numCombinations
                                for n=1:numElectrodes
                                    spikeElectrode=spikeElectrodes(n);
        
                                    disp(['Working on electrode: ' num2str(spikeElectrode)]);
                                    tmp = load(fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(spikeElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality{moda} '_' tuning 'tuning_Allmethods_for_distance' electrodeDistanceVal{dist} '.mat']));
                                    blah=cell2mat(tmp.staVals);
                                
                                                                 
                                   for ee=1:size(blah,1)
                                        for tt=1:size(blah,4)                       
                                            staValsMat(ee,1,1,tt,:)=squeeze(blah(ee,:,1,tt));
                                        end
                                   end
        
                                   for ee=1:size(tmp.ffc,1) % number of LFP elec
                                    for tt=1:size(tmp.ffc,4) % number of analysis period                           
                        
                                        lfp_ffc(n,ee,1,1,tt,:)=tmp.ffc(ee,1,1,tt,:);
                                        lfp_sfc(n,ee,1,1,tt,:)=tmp.sfc(ee,1,1,tt,:);
                                        lfp_staValsMat(n,ee,1,1,tt,:)=staValsMat(ee,1,1,tt,:);
                                    end
                                   end
                                end    
        
                                ffcData{dist}(1,:)=squeeze(mean(lfp_ffc(:,:,1,1,1,:),[1 2])); %average across all spike and lfp neighnboring electrodes in baseline
                                ffcData{dist}(2,:)=squeeze(mean(lfp_ffc(:,:,1,1,2,:),[1 2])); %average across all spike and lfp neighnboring electrodes in stim
    
                                sfcData{dist}(1,:)=squeeze(mean(lfp_sfc(:,:,1,1,1,:),[1 2])); %average across all spike and lfp neighnboring electrodes in baseline
                                sfcData{dist}(2,:)=squeeze(mean(lfp_sfc(:,:,1,1,2,:),[1 2])); %average across all spike and lfp neighnboring electrodes in stim
    
                                staData{dist}(1,:)=squeeze(mean(lfp_staValsMat(:,:,1,1,1,:),[1 2])); %average across all spike and lfp neighnboring electrodes in baseline
                                staData{dist}(2,:)=squeeze(mean(lfp_staValsMat(:,:,1,1,2,:),[1 2])); %average across all spike and lfp neighnboring electrodes in stim
                            end
                end
            end                


            %plotting
            ffcFreq = tmp.ffcFreq;
            sfcFreq = tmp.sfcFreq;
            xsSTA = 1000*tmp.xsSTA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for o=1:numOrientations
                for i=1:numCombinations
                    
                    plot(plotHandlesFFC(o),ffcFreq,ffcData{i}(o,:),'color',distanceColorNames(i));
                    hold(plotHandlesFFC(o),'on');
            
                    plot(plotHandlesSFC(o),sfcFreq,sfcData{i}(o,:),'color',distanceColorNames(i));
                    hold(plotHandlesSFC(o),'on');
            
                    plot(plotHandlesSTA(o),xsSTA,staData{i}(o,:),'color',distanceColorNames(i));
                    hold(plotHandlesSTA(o),'on');
        
                    
                end
            
                set(plotHandlesFFC(o),'XTickLabel',[]);
                if o==1
                    ylabel(plotHandlesFFC(o),'FFC');
                    ylabel(plotHandlesSFC(o),'SFC');
                    ylabel(plotHandlesSTA(o),'STA');
                    title(plotHandlesFFC(o),[labelNames{o}]);
                    
                else
                    set(plotHandlesFFC(o),'YTickLabel',[]);
                    set(plotHandlesSFC(o),'YTickLabel',[]);
                    set(plotHandlesSTA(o),'YTickLabel',[]);
        
                    title(plotHandlesFFC(o),[labelNames{o}]);
                    
                end
            
                xlabel(plotHandlesSFC(o),'Frequency (Hz)');
                xlabel(plotHandlesSTA(o),'Time (ms)');
        
                axis(plotHandlesFFC(o),[0 ffcFreq(end) 0 1]);
                axis(plotHandlesSFC(o),[0 ffcFreq(end) 0 1]);
                axis(plotHandlesSTA(o),[xsSTA(1) xsSTA(end) -80 40]);
                
        
            end
            

        % plotting either PAC values or corresponding p-values
        elseif (strcmp(metrics,'pac') || strcmp(metrics,'pval'))
                       

            labelNames={'Morlet-bsl', 'Morlet-stim','FiltFilt-bsl', 'FiltFilt-stim','MP-bsl','MP-stim',};
            numOrientations=6; 
           

            plotHandlesPAC1 = getPlotHandles(1,numOrientations,[0.1 0.75 0.8 0.11],0); linkaxes(plotHandlesPAC1);
            plotHandlesPAC2 = getPlotHandles(1,numOrientations,[0.1 0.60 0.8 0.11],0); linkaxes(plotHandlesPAC2);
            plotHandlesPAC3 = getPlotHandles(1,numOrientations,[0.1 0.45 0.8 0.11],0); linkaxes(plotHandlesPAC3);
            plotHandlesPAC4 = getPlotHandles(1,numOrientations,[0.1 0.30 0.8 0.11],0); linkaxes(plotHandlesPAC4);
            plotHandlesPAC5 = getPlotHandles(1,numOrientations,[0.1 0.15 0.8 0.11],0); linkaxes(plotHandlesPAC5);
            

            switch modality

                case 'ECoG'
            
                    if strcmp(monkeyName,'alpaH')
            
                        ecogElectrodes=[82 85 86 88 89];
            
                    elseif strcmp(monkeyName,'kesariH')
            
                        ecogElectrodes=[85 86 88 89];
                    end
                    lfpElectrodes=ecogElectrodes;
                    
   
                    numElectrodes=length(ecogElectrodes);
                    
                    if strcmp(electrodeDistanceVal,'400')
                        disp('file not found');   

                    else
                        for n=1:numElectrodes
                            ecogElectrode=ecogElectrodes(n);

                            disp(['Working on electrode: ' num2str(ecogElectrode)]);
                            tmp = load(fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(ecogElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_' tuning 'tuning_Allmethods_for_distance' electrodeDistanceVal '.mat']));


                           for ee=1:size(tmp.pacmat,1) % number of LFP elec
                            for tt=1:size(tmp.pacmat,4) % number of analysis period    
                                for fi=1:size(tmp.pacmat,5) % number of filter method
                                    for meth=1:size(tmp.pacmat,6) % number of filter method
                
                                        pacval=squeeze(tmp.pacmat(ee,1,1,tt,fi,meth,:,:));
                                        pval=squeeze(tmp.pval(ee,1,1,tt,fi,meth,:,:));
                                        pacval(isnan(pacval))=0;
                                        if ~isreal(pacval)
                                            pacval=abs(pacval);
                                        end
                                        clear corrected_pval survived_pacval
                                        %performing p value correction
                                        %either FDR or Bonferrroni or
                                        %remain uncorrected
                                        if strcmp(correction,'fdr')
                                            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval,pThresh,'pdep','yes');
                                            corrected_pval=adj_p.*h;
                                            survived_pacval=pacval.*h;
                                        elseif strcmp(correction,'bonf')
                                            pThreshNew=pThresh/(size(pacval,1)*size(pacval,2));
                                            idx=find(pval<pThreshNew);
                                            corrected_pval=pval*0;
                                            corrected_pval(idx)=pval(idx);
                                            survived_pacval=pacval*0;
                                            survived_pacval(idx)=pacval(idx);
                                        elseif strcmp(correction,'none')
                                            corrected_pval=pval;
                                            survived_pacval=pacval;
                                        end
                                        pacmatAll(n,ee,1,1,tt,fi,meth,:,:)=survived_pacval;
                                        pvalAll(n,ee,1,1,tt,fi,meth,:,:)=corrected_pval;
                                    end
                                end
                            end
                           end
                        end                        
                    end
                                    
            
                case 'LFP'
            
                    lfpElectrodes = goodElectrodes(goodElectrodes<=81);
                    numElectrodes=length(spikeElectrodes);
                    for n=1:numElectrodes
                        spikeElectrode=spikeElectrodes(n);
    
                        disp(['Working on electrode: ' num2str(spikeElectrode)]);
                        tmp = load(fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(spikeElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_' tuning 'tuning_Allmethods_for_distance' electrodeDistanceVal '.mat']));
    
    
                       for ee=1:size(tmp.pacmat,1) % number of LFP elec
                        for tt=1:size(tmp.pacmat,4) % number of analysis period    
                            for fi=1:size(tmp.pacmat,5) % number of filter method
                                for meth=1:size(tmp.pacmat,6) % number of filter method
            
                                    pacval=squeeze(tmp.pacmat(ee,1,1,tt,fi,meth,:,:));
                                    pval=squeeze(tmp.pval(ee,1,1,tt,fi,meth,:,:));
                                    pacval(isnan(pacval))=0;
                                    if ~isreal(pacval)
                                        pacval=abs(pacval);
                                    end
                                    clear corrected_pval survived_pacval
                                    %performing p value correction
                                    %either FDR or Bonferrroni or
                                    %remain uncorrected
                                    if strcmp(correction,'fdr')
                                        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval,pThresh,'pdep','yes');
                                        corrected_pval=adj_p.*h;
                                        survived_pacval=pacval.*h;
                                    elseif strcmp(correction,'bonf')
                                        pThreshNew=pThresh/(size(pacval,1)*size(pacval,2));
                                        idx=find(pval<pThreshNew);
                                        corrected_pval=pval*0;
                                        corrected_pval(idx)=pval(idx);
                                        survived_pacval=pacval*0;
                                        survived_pacval(idx)=pacval(idx);
                                    elseif strcmp(correction,'none')
                                        corrected_pval=pval;
                                        survived_pacval=pacval;
                                        
                                    end
                                    pacmatAll(n,ee,1,1,tt,fi,meth,:,:)=survived_pacval;
                                    pvalAll(n,ee,1,1,tt,fi,meth,:,:)=corrected_pval;
                                end
                            end
                        end
                       end
                    end                        
            end

            avgPacData(1,:,:,:)=squeeze(mean(pacmatAll(:,:,1,1,1,1,:,:,:),[ 1 2])); % bsl Filt1
            avgPacData(2,:,:,:)=squeeze(mean(pacmatAll(:,:,1,1,2,1,:,:,:),[ 1 2])); % stim Filt1
            avgPacData(3,:,:,:)=squeeze(mean(pacmatAll(:,:,1,1,1,2,:,:,:),[ 1 2])); % bsl Filt2
            avgPacData(4,:,:,:)=squeeze(mean(pacmatAll(:,:,1,1,2,2,:,:,:),[ 1 2])); % stim Filt2
            avgPacData(5,:,:,:)=squeeze(mean(pacmatAll(:,:,1,1,1,3,:,:,:),[ 1 2])); % bsl Filt3
            avgPacData(6,:,:,:)=squeeze(mean(pacmatAll(:,:,1,1,2,3,:,:,:),[ 1 2])); % stim Filt3

            avgPvalData(1,:,:,:)=squeeze(mean(pvalAll(:,:,1,1,1,1,:,:,:),[ 1 2])); % bsl Filt1
            avgPvalData(2,:,:,:)=squeeze(mean(pvalAll(:,:,1,1,2,1,:,:,:),[ 1 2])); % stim Filt1
            avgPvalData(3,:,:,:)=squeeze(mean(pvalAll(:,:,1,1,1,2,:,:,:),[ 1 2])); % bsl Filt2
            avgPvalData(4,:,:,:)=squeeze(mean(pvalAll(:,:,1,1,2,2,:,:,:),[ 1 2])); % stim Filt2
            avgPvalData(5,:,:,:)=squeeze(mean(pvalAll(:,:,1,1,1,3,:,:,:),[ 1 2])); % bsl Filt3
            avgPvalData(6,:,:,:)=squeeze(mean(pvalAll(:,:,1,1,2,3,:,:,:),[ 1 2])); % stim Filt3                                      


            %plotting
            PhaseFreqVector=tmp.freqvec_ph;
            AmpFreqVector=tmp.freqvec_amp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(metrics,'pac')
                dataPlot=avgPacData;
            elseif strcmp(metrics,'pval')
                dataPlot=-log(avgPvalData);
                dataPlot(isinf(dataPlot))=0;
            end

            for o=1:numOrientations
                    
                contourf(plotHandlesPAC1(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,1,:,:)) ,30,'lines','none');  
                hold(plotHandlesPAC1(o),'on');

                contourf(plotHandlesPAC2(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,2,:,:)),30,'lines','none');  
                hold(plotHandlesPAC2(o),'on');

                contourf(plotHandlesPAC3(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,3,:,:)),30,'lines','none');  
                hold(plotHandlesPAC3(o),'on');

                contourf(plotHandlesPAC4(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,4,:,:)),30,'lines','none');  
                hold(plotHandlesPAC4(o),'on');

                contourf(plotHandlesPAC5(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,5,:,:)),30,'lines','none');  
                hold(plotHandlesPAC5(o),'on');                  
                    
                
            
                set(plotHandlesPAC1(o),'XTickLabel',[]);
                set(plotHandlesPAC2(o),'XTickLabel',[]);
                set(plotHandlesPAC3(o),'XTickLabel',[]);
                set(plotHandlesPAC4(o),'XTickLabel',[]);
                if o==1
                    ylabel(plotHandlesPAC5(o),sprintf('Amplitude Frequency (Hz)'));
                    title(plotHandlesPAC1(o),[labelNames{o}],'FontWeight','bold','FontSize',12);
                    
                else
                    set(plotHandlesPAC1(o),'YTickLabel',[]);
                    set(plotHandlesPAC2(o),'YTickLabel',[]);
                    set(plotHandlesPAC3(o),'YTickLabel',[]);
                    set(plotHandlesPAC4(o),'YTickLabel',[]);
                    set(plotHandlesPAC5(o),'YTickLabel',[]);
                    title(plotHandlesPAC1(o),[labelNames{o}],'FontWeight','bold','FontSize',12);
                            
                end


             if rem(o,2)==1
                valmin1=min(squeeze(dataPlot(o:o+1,1,:,:)),[],'all'); valmax1=max(squeeze(dataPlot(o:o+1,1,:,:)),[],'all');
                clim(plotHandlesPAC1(o),[ valmin1 valmax1]);

                valmin2=min(squeeze(dataPlot(o:o+1,2,:,:)),[],'all'); valmax2=max(squeeze(dataPlot(o:o+1,2,:,:)),[],'all');
                clim(plotHandlesPAC2(o),[ valmin2 valmax2]);

                valmin3=min(squeeze(dataPlot(o:o+1,3,:,:)),[],'all'); valmax3=max(squeeze(dataPlot(o:o+1,3,:,:)),[],'all');
                clim(plotHandlesPAC3(o),[ valmin3 valmax3]);

                valmin4=min(squeeze(dataPlot(o:o+1,4,:,:)),[],'all'); valmax4=max(squeeze(dataPlot(o:o+1,4,:,:)),[],'all');
                clim(plotHandlesPAC4(o),[ valmin4 valmax4]);

                valmin5=min(squeeze(dataPlot(o:o+1,5,:,:)),[],'all'); valmax5=max(squeeze(dataPlot(o:o+1,5,:,:)),[],'all');
                clim(plotHandlesPAC5(o),[ valmin5 valmax5]);

            elseif rem(o,2)==0
                valmin1=min(squeeze(dataPlot(o-1:o,1,:,:)),[],'all'); 
                valmax1=max(squeeze(dataPlot(o-1:o,1,:,:)),[],'all');
                clim(plotHandlesPAC1(o),[ valmin1 valmax1]);
%                 colormap(plotHandlesPAC1(o),"autumn");            
                cb=colorbar(plotHandlesPAC1(o),'south');
                cb.Position = cb.Position + [-0.07 -0.05 0 0];
                cb.Ticks = [valmin1 valmax1];
                x1=xlabel(cb,'MVLMI','FontWeight','bold','FontSize',8,'Rotation',0);
                x1.Position(2) = min(ylim(cb)) ;x1.VerticalAlignment="middle";
                x1.Position(2) = min(ylim(cb))+1.8;


                valmin2=min(squeeze(dataPlot(o-1:o,2,:,:)),[],'all'); 
                valmax2=max(squeeze(dataPlot(o-1:o,2,:,:)),[],'all');
                clim(plotHandlesPAC2(o),[ valmin2 valmax2]);
%                 colormap(plotHandlesPAC1(o),"autumn");            
                cb=colorbar(plotHandlesPAC2(o),'south');
                cb.Position = cb.Position + [-0.07 -0.05 0 0];
                cb.Ticks = [valmin2 valmax2];
                x1=xlabel(cb,'PLV','FontWeight','bold','FontSize',8,'Rotation',0);
                x1.Position(2) = min(ylim(cb)) ;x1.VerticalAlignment="middle";
                x1.Position(2) = min(ylim(cb))+1.8;

                valmin3=min(squeeze(dataPlot(o-1:o,3,:,:)),[],'all'); 
                valmax3=max(squeeze(dataPlot(o-1:o,3,:,:)),[],'all');
                clim(plotHandlesPAC3(o),[ valmin3 valmax3]);
%                 colormap(plotHandlesPAC1(o),"autumn");            
                cb=colorbar(plotHandlesPAC3(o),'south');
                cb.Position = cb.Position + [-0.07 -0.05 0 0];
                cb.Ticks = [valmin3 valmax3];
                x1=xlabel(cb,'GLM','FontWeight','bold','FontSize',8,'Rotation',0);
                x1.Position(2) = min(ylim(cb)) ;x1.VerticalAlignment="middle";
                x1.Position(2) = min(ylim(cb))+1.8;

                valmin4=min(squeeze(dataPlot(o-1:o,4,:,:)),[],'all'); 
                valmax4=max(squeeze(dataPlot(o-1:o,4,:,:)),[],'all');
                clim(plotHandlesPAC4(o),[ valmin4 valmax4]);
%                 colormap(plotHandlesPAC1(o),"autumn");            
                cb=colorbar(plotHandlesPAC4(o),'south');
                cb.Position = cb.Position + [-0.07 -0.05 0 0];
                cb.Ticks = [valmin4 valmax4];
                x1=xlabel(cb,'KLMI','FontWeight','bold','FontSize',8,'Rotation',0);
                x1.Position(2) = min(ylim(cb)) ;x1.VerticalAlignment="middle";
                x1.Position(2) = min(ylim(cb))+1.8;

                valmin5=min(squeeze(dataPlot(o-1:o,5,:,:)),[],'all'); 
                valmax5=max(squeeze(dataPlot(o-1:o,5,:,:)),[],'all');
                clim(plotHandlesPAC5(o),[ valmin5 valmax5]);
%                 colormap(plotHandlesPAC1(o),"autumn");            
                cb=colorbar(plotHandlesPAC5(o),'south');
                cb.Position = cb.Position + [-0.07 -0.1 0 0];
                cb.Ticks = [valmin5 valmax5];
                x1=xlabel(cb,'n-MVLMI','FontWeight','bold','FontSize',8,'Rotation',0);
                x1.Position(2) = min(ylim(cb)) ;x1.VerticalAlignment="middle";
                x1.Position(2) = min(ylim(cb))+1.8;
            end
        
            end

            xlabel(plotHandlesPAC5(1),'Phase Frequency from LFP data (Hz)');      
           
        end      

    case 'multiple' % for SF, Ori and size tuning with multiple conditions

        numElectrodes = length(spikeElectrodes);
        if strcmp(tuning,'Ori')
            numOrientations = 9;            
        elseif strcmp(tuning,'SF')
            numOrientations=6;
        elseif strcmp(tuning,'Size')
            numOrientations=7;
        end
        for nfil=1:length(filterMethod)
        
            figure(nfil) %seperate figure for each filtering methods

            plotHandlesFFC = getPlotHandles(1,numOrientations,[0.05 0.85 0.9 0.125],0.01); linkaxes(plotHandlesFFC);
            plotHandlesSFC = getPlotHandles(1,numOrientations,[0.05 0.70 0.9 0.125],0.01); linkaxes(plotHandlesSFC);
            plotHandlesSTA = getPlotHandles(1,numOrientations,[0.05 0.55 0.9 0.1],0.01); linkaxes(plotHandlesSTA);
            plotHandlesPAC1 = getPlotHandles(1,numOrientations,[0.05 0.40 0.9 0.06],0.01); linkaxes(plotHandlesPAC1);
            plotHandlesPAC2 = getPlotHandles(1,numOrientations,[0.05 0.32 0.9 0.06],0.01); linkaxes(plotHandlesPAC2);
            plotHandlesPAC3 = getPlotHandles(1,numOrientations,[0.05 0.24 0.9 0.06],0.01); linkaxes(plotHandlesPAC3);
            plotHandlesPAC4 = getPlotHandles(1,numOrientations,[0.05 0.16 0.9 0.06],0.01); linkaxes(plotHandlesPAC4);
            plotHandlesPAC5 = getPlotHandles(1,numOrientations,[0.05 0.08 0.9 0.06],0.01); linkaxes(plotHandlesPAC5);
    
    
    
            for n=1:numElectrodes
                spikeElectrode = spikeElectrodes(n);
            
                disp(['Working on electrode: ' num2str(spikeElectrode)]);
                tmp = load(fullfile(folderSave,[monkeyName expDate protocolName 'elec' num2str(spikeElectrode) '_removeMean' num2str(removeEvokedResponse) 'Tapers' num2str(tapers(1)) '_' num2str(tapers(2)) '_' modality '_' tuning 'tuning_Allmethods_for_distance' electrodeDistanceVal '.mat']));
            
                blah=cell2mat(tmp.staVals);
            
                tmp.pacmat(isnan(tmp.pacmat))=0;
            
            %     clear staValsMat
                if strcmp(tuning,'Ori')
                   for ee=1:size(blah,1)
                    for cc=1:size(blah,3)
                        for tt=1:size(blah,4)                       
                            staValsMat(ee,1,cc,tt,:)=squeeze(blah(ee,:,cc,tt));
                        end
                     end
                   end                  
           
            
                    for ee=1:size(tmp.ffc,1)
                        for tt=1:size(tmp.ffc,4)
                            for nOri=1:numOrientations
                                % Since 'all' was excluded during computation, while plotting generating 'all' condition by averaging across SF or Ori pac values            
                                if nOri==9
                                    tmp_ffc(n,ee,1,9,tt,:)=squeeze(mean(tmp.ffc(ee,1,1:8,tt,:),3));
                                    tmp_sfc(n,ee,1,9,tt,:)=squeeze(mean(tmp.sfc(ee,1,1:8,tt,:),3));
                                    tmp_staValsMat(n,ee,1,9,tt,:)=squeeze(mean(staValsMat(ee,1,1:8,tt,:),3));
                                    tmp_pacmat(n,ee,1,9,tt,:,:,:,:)=squeeze(mean(tmp.pacmat(ee,1,1:8,tt,:,:,:,:),3));
                                    tmp_pval(n,ee,1,9,tt,:,:,:,:)=squeeze(mean(tmp.pval(ee,1,1:8,tt,:,:,:,:),3));
                                else
                                    tmp_ffc(n,ee,1,nOri,tt,:)=tmp.ffc(ee,1,nOri,tt,:);
                                    tmp_sfc(n,ee,1,nOri,tt,:)=tmp.sfc(ee,1,nOri,tt,:);
                                    tmp_staValsMat(n,ee,1,nOri,tt,:)=staValsMat(ee,1,nOri,tt,:);
                                    tmp_pacmat(n,ee,1,nOri,tt,:,:,:,:)=tmp.pacmat(ee,1,nOri,tt,:,:,:,:);
                                    tmp_pval(n,ee,1,nOri,tt,:,:,:,:)=tmp.pval(ee,1,nOri,tt,:,:,:,:);
                                end
                            end
                        end
                    end
    
                elseif strcmp(tuning,'SF')
                    len=size(tmp.staVals{1},2);
                    for ee=1:size(tmp.staVals,1)
                            for cc=1:size(tmp.staVals,2)
                                for tt=1:size(tmp.staVals,4) 
                                    if cc==1
                                        staValsMat(ee,cc,1,tt,:)=squeeze(blah(ee,cc:cc*len,1,tt));
                                    else
                                        staValsMat(ee,cc,1,tt,:)=squeeze(blah(ee,(1+len*(cc-1)):(cc*len),1,tt));
                                    end
                                end
                            end
                    end
                    
                    for ee=1:size(tmp.ffc,1)
                        for tt=1:size(tmp.ffc,4)
                            for nOri=1:numOrientations
                               % Since 'all' was excluded during computation, while plotting generating 'all' condition by averaging across SF or Ori pac values            
                                if nOri==6
                                    tmp_ffc(n,ee,6,1,tt,:)=squeeze(mean(tmp.ffc(ee,1:5,1,tt,:),2));
                                    tmp_sfc(n,ee,6,1,tt,:)=squeeze(mean(tmp.sfc(ee,1:5,1,tt,:),2));
                                    tmp_staValsMat(n,ee,6,1,tt,:)=squeeze(mean(staValsMat(ee,1:5,1,tt,:),2));
                                    tmp_pacmat(n,ee,6,1,tt,:,:,:,:)=squeeze(mean(tmp.pacmat(ee,1:5,1,tt,:,:,:,:),2));
                                    tmp_pval(n,ee,6,1,tt,:,:,:,:)=squeeze(mean(tmp.pval(ee,1:5,1,tt,:,:,:,:),2));

                                else
                                    tmp_ffc(n,ee,nOri,1,tt,:)=tmp.ffc(ee,nOri,1,tt,:);
                                    tmp_sfc(n,ee,nOri,1,tt,:)=tmp.sfc(ee,nOri,1,tt,:);
                                    tmp_staValsMat(n,ee,nOri,1,tt,:)=staValsMat(ee,nOri,1,tt,:);
                                    tmp_pacmat(n,ee,nOri,1,tt,:,:,:,:)=squeeze(tmp.pacmat(ee,nOri,1,tt,:,:,:,:));
                                    tmp_pval(n,ee,nOri,1,tt,:,:,:,:)=squeeze(tmp.pval(ee,nOri,1,tt,:,:,:,:));

                                end
                            end     
                        end
                    end
                elseif strcmp(tuning,'Size')
                    len=size(tmp.staVals{1},2);
                    for ee=1:size(tmp.staVals,1)
                            for cc=1:size(tmp.staVals,2)
                                for tt=1:size(tmp.staVals,4) 
                                    if cc==1
                                        staValsMat(ee,cc,1,tt,:)=squeeze(blah(ee,cc:cc*len,1,tt));
                                    else
                                        staValsMat(ee,cc,1,tt,:)=squeeze(blah(ee,(1+len*(cc-1)):(cc*len),1,tt));
                                    end
                                end
                            end
                    end
                    
                    for ee=1:size(tmp.ffc,1)
                        for tt=1:size(tmp.ffc,4)
                            for nOri=1:numOrientations
                             % Since 'all' was excluded during computation, while plotting generating 'all' condition by averaging across Size pac values            

                                if nOri==7
                                    tmp_ffc(n,ee,7,1,tt,:)=squeeze(mean(tmp.ffc(ee,1:6,1,tt,:),2));
                                    tmp_sfc(n,ee,7,1,tt,:)=squeeze(mean(tmp.sfc(ee,1:6,1,tt,:),2));
                                    tmp_staValsMat(n,ee,7,1,tt,:)=squeeze(mean(staValsMat(ee,1:6,1,tt,:),2));
                                    tmp_pacmat(n,ee,7,1,tt,:,:,:,:)=squeeze(mean(tmp.pacmat(ee,1:6,1,tt,:,:,:,:),2));
                                    tmp_pval(n,ee,7,1,tt,:,:,:,:)=squeeze(mean(tmp.pval(ee,1:6,1,tt,:,:,:,:),2));

                                else
                                    tmp_ffc(n,ee,nOri,1,tt,:)=tmp.ffc(ee,nOri,1,tt,:);
                                    tmp_sfc(n,ee,nOri,1,tt,:)=tmp.sfc(ee,nOri,1,tt,:);
                                    tmp_staValsMat(n,ee,nOri,1,tt,:)=staValsMat(ee,nOri,1,tt,:);
                                    tmp_pacmat(n,ee,nOri,1,tt,:,:,:,:)=squeeze(tmp.pacmat(ee,nOri,1,tt,:,:,:,:));
                                    tmp_pval(n,ee,nOri,1,tt,:,:,:,:)=squeeze(tmp.pval(ee,nOri,1,tt,:,:,:,:));

                                end
                            end     
                        end
                    end



            
                end
            end
    
    % performing correction and getting only corrected pacval and pval              
    
            for n=1:size(tmp_pacmat,1) %number of spikeElectrode
              for ee=1:size(tmp_pacmat,2) % number of LFP elec
                 for nSF=1:size(tmp_pacmat,3) % number of SF
                     for nOri=1:size(tmp_pacmat,4) % number of Ori
                         for tt=1:size(tmp_pacmat,5) % number of analysis period   
                             for fi=1:size(tmp_pacmat,6) % number of filter method
                                 for meth=1:size(tmp_pacmat,7) % number of PAC method
    
                                    pacval=squeeze(tmp_pacmat(n,ee,nSF,nOri,tt,fi,meth,:,:));
                                    pval=squeeze(tmp_pval(n,ee,nSF,nOri,tt,fi,meth,:,:));
                                    pacval(isnan(pacval))=0;
                                    if ~isreal(pacval)
                                        pacval=abs(pacval);
                                    end
                                    clear corrected_pval survived_pacval
                                    if strcmp(correction,'fdr')
                                        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval,pThresh,'pdep','yes');
                                        corrected_pval=adj_p.*h;
                                        survived_pacval=pacval.*h;
                                    elseif strcmp(correction,'bonf')
                                        pThreshNew=pThresh/(size(pacval,1)*size(pacval,2));
                                        idx=find(pval<pThreshNew);
                                        corrected_pval=pval*0;
                                        corrected_pval(idx)=pval(idx);
                                        survived_pacval=pacval*0;
                                        survived_pacval(idx)=pacval(idx);
                                    elseif strcmp(correction,'none')
                                        corrected_pval=pval;
                                        survived_pacval=pacval;
                                        
                                    end
                                    pacmatAll(n,ee,nSF,nOri,tt,fi,meth,:,:)=survived_pacval;
                                    pvalAll(n,ee,nSF,nOri,tt,fi,meth,:,:)=corrected_pval;
                                 end
                             end
                         end
                     end
                 end
              end
            end            
    
        if strcmp(tuning,'Ori')   
            for j=1:numCombinations
                ffcData{j} = squeeze(mean(tmp_ffc(:,:,1,:,analysisPeriodNum,:),[1 2]));
                sfcData{j} =squeeze(mean(tmp_sfc(:,:,1,:,analysisPeriodNum,:),[1 2]));
                staData{j} = squeeze(mean(tmp_staValsMat(:,:,1,:,analysisPeriodNum,:),[1 2]));
                avgPacData=squeeze(mean(pacmatAll(:,:,1,:,analysisPeriodNum,nfil,:,:,:),[1 2]));
                avgPvalData=squeeze(mean(pvalAll(:,:,1,:,analysisPeriodNum,nfil,:,:,:),[1 2]));
    
            end
        elseif strcmp(tuning,'SF') || strcmp(tuning,'Size')
            for j=1:numCombinations
                ffcData{j} = squeeze(mean(tmp_ffc(:,:,:,1,analysisPeriodNum,:),[1 2]));
                sfcData{j} = squeeze(mean(tmp_sfc(:,:,:,1,analysisPeriodNum,:),[1 2]));
                staData{j} = squeeze(mean(tmp_staValsMat(:,:,:,1,analysisPeriodNum,:),[1 2]));
                avgPacData=squeeze(mean(pacmatAll(:,:,:,1,analysisPeriodNum,nfil,:,:,:),[1 2]));
                avgPvalData=squeeze(mean(pvalAll(:,:,:,1,analysisPeriodNum,nfil,:,:,:),[1 2]));
    
            end
    
        end
    
    ffcFreq = tmp.ffcFreq;
    sfcFreq = tmp.sfcFreq;
    xsSTA = 1000*tmp.xsSTA;
    PhaseFreqVector=tmp.freqvec_ph;
    AmpFreqVector=tmp.freqvec_amp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(metrics,'pac')
        dataPlot=avgPacData;
    elseif strcmp(metrics,'pval')
        dataPlot=-log(avgPvalData);
        dataPlot(isinf(dataPlot))=0;
    end 

    if strcmp(tuning,'Ori')
        oriNames={'0' '22' '45' '67' '90' '112' '135' '157' 'all'};
    
        for o=1:numOrientations
            for i=1:numCombinations
                
                plot(plotHandlesFFC(o),ffcFreq,ffcData{i}(o,:),'color',distanceColorNames(i));
                hold(plotHandlesFFC(o),'on');
    
        
                plot(plotHandlesSFC(o),sfcFreq,sfcData{i}(o,:),'color',distanceColorNames(i));
                hold(plotHandlesSFC(o),'on');
        
                plot(plotHandlesSTA(o),xsSTA,staData{i}(o,:),'color',distanceColorNames(i));
                hold(plotHandlesSTA(o),'on');
                    
                contourf(plotHandlesPAC1(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,1,:,:)),30,'lines','none');  
                hold(plotHandlesPAC1(o),'on');
    
                contourf(plotHandlesPAC2(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,2,:,:)),30,'lines','none');  
                hold(plotHandlesPAC2(o),'on');
    
                contourf(plotHandlesPAC3(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,3,:,:)),30,'lines','none');  
                hold(plotHandlesPAC3(o),'on');
    
                contourf(plotHandlesPAC4(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,4,:,:)),30,'lines','none');  
                hold(plotHandlesPAC4(o),'on');
    
                contourf(plotHandlesPAC5(o),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(o,5,:,:)),30,'lines','none');  
                hold(plotHandlesPAC5(o),'on');
                
            end
        
            set(plotHandlesFFC(o),'XTickLabel',[]);
            if o==1
                ylabel(plotHandlesFFC(o),'FFC');
                ylabel(plotHandlesSFC(o),'SFC');
                ylabel(plotHandlesSTA(o),'STA');
                ylabel(plotHandlesPAC5(o),sprintf('Amplitude Frequency (Hz)'));
                title(plotHandlesFFC(o),['SF=4 Ori=' oriNames{o}]);
                
            else
                set(plotHandlesFFC(o),'YTickLabel',[]);
                set(plotHandlesSFC(o),'YTickLabel',[]);
                set(plotHandlesSTA(o),'YTickLabel',[]);
                set(plotHandlesPAC1(o),'YTickLabel',[]);
                set(plotHandlesPAC2(o),'YTickLabel',[]);
                set(plotHandlesPAC3(o),'YTickLabel',[]);
                set(plotHandlesPAC4(o),'YTickLabel',[]);
                set(plotHandlesPAC5(o),'YTickLabel',[]);
                title(plotHandlesFFC(o),['SF=4 Ori=' oriNames{o}]);
    
    
                
            end
        
            xlabel(plotHandlesSFC(o),'Frequency (Hz)');
            xlabel(plotHandlesSTA(o),'Time (ms)');
            xlabel(plotHandlesPAC5(1),'Phase Frequency from LFP data (Hz)');
    
            axis(plotHandlesFFC(o),[0 ffcFreq(end) 0 1]);
            axis(plotHandlesSFC(o),[0 ffcFreq(end) 0 1]);
            axis(plotHandlesSTA(o),[xsSTA(1) xsSTA(end) -80 40]);
            
            clim(plotHandlesPAC1(o),[min(dataPlot(:,1,:,:),[],'all') max(dataPlot(:,1,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC1(9));
            cb.Ticks = [min(dataPlot(:,1,:,:),[],'all') max(dataPlot(:,1,:,:),[],'all')];
            ylabel(cb,'MVLMI','FontWeight','bold','FontSize',8,'Rotation',270)
    
            clim(plotHandlesPAC2(o),[min(dataPlot(o,2,:,:),[],'all') max(dataPlot(:,2,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC2(9));
            cb.Ticks = [min(dataPlot(:,2,:,:),[],'all') max(dataPlot(:,2,:,:),[],'all')];
            ylabel(cb,'PLV','FontWeight','bold','FontSize',8,'Rotation',270)
    
            clim(plotHandlesPAC3(o),[min(dataPlot(:,3,:,:),[],'all') max(dataPlot(:,3,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC3(9));
            cb.Ticks = [min(dataPlot(:,3,:,:),[],'all') max(dataPlot(:,3,:,:),[],'all')];
            ylabel(cb,'GLM','FontWeight','bold','FontSize',8,'Rotation',270)
    
    
            clim(plotHandlesPAC4(o),[min(dataPlot(:,4,:,:),[],'all') max(dataPlot(:,4,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC4(9));
            cb.Ticks = [min(dataPlot(:,4,:,:),[],'all') max(dataPlot(o,4,:,:),[],'all')];
            ylabel(cb,'KLMI','FontWeight','bold','FontSize',8,'Rotation',270)
    
    
            clim(plotHandlesPAC5(o),[min(dataPlot(:,5,:,:),[],'all') max(dataPlot(:,5,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC5(9));
            cb.Ticks = [min(dataPlot(o,5,:,:),[],'all') max(dataPlot(o,5,:,:),[],'all')];
            ylabel(cb,'n-MVL','FontWeight','bold','FontSize',8,'Rotation',270)
    
    
    
        end
    elseif strcmp(tuning,'SF')
            sfNames={'0.5', '1', '2', '4', '8', 'all'};
            for f=1:numOrientations
            for i=1:numCombinations
                
                plot(plotHandlesFFC(f),ffcFreq,ffcData{i}(f,:),'color',distanceColorNames(i));
                hold(plotHandlesFFC(f),'on');
        
                plot(plotHandlesSFC(f),sfcFreq,sfcData{i}(f,:),'color',distanceColorNames(i));
                hold(plotHandlesSFC(f),'on');
        
                plot(plotHandlesSTA(f),xsSTA,staData{i}(f,:),'color',distanceColorNames(i));
                hold(plotHandlesSTA(f),'on');
                
                contourf(plotHandlesPAC1(f),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(f,1,:,:)),30,'lines','none');  
                hold(plotHandlesPAC1(f),'on');

                contourf(plotHandlesPAC2(f),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(f,2,:,:)),30,'lines','none');  
                hold(plotHandlesPAC2(f),'on');

                contourf(plotHandlesPAC3(f),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(f,3,:,:)),30,'lines','none');  
                hold(plotHandlesPAC3(f),'on');

                contourf(plotHandlesPAC4(f),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(f,4,:,:)),30,'lines','none');  
                hold(plotHandlesPAC4(f),'on');

                contourf(plotHandlesPAC5(f),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(f,5,:,:)),30,'lines','none');  
                hold(plotHandlesPAC5(f),'on');
    
            end
        
            set(plotHandlesFFC(f),'XTickLabel',[]);
            if f==1
                ylabel(plotHandlesFFC(f),'FFC');
                ylabel(plotHandlesSFC(f),'SFC');
                ylabel(plotHandlesSTA(f),'STA');           
                
                ylabel(plotHandlesPAC5(f),sprintf('Amplitude Frequency (Hz)'));
                title(plotHandlesFFC(f),['SF=' sfNames{f} 'Ori=90' ]);
    
            else
                set(plotHandlesFFC(f),'YTickLabel',[]);
                set(plotHandlesSFC(f),'YTickLabel',[]);
                set(plotHandlesSTA(f),'YTickLabel',[]);
                set(plotHandlesPAC1(f),'YTickLabel',[]);
                set(plotHandlesPAC2(f),'YTickLabel',[]);
                set(plotHandlesPAC3(f),'YTickLabel',[]);
                set(plotHandlesPAC4(f),'YTickLabel',[]);
                set(plotHandlesPAC5(f),'YTickLabel',[]);
                title(plotHandlesFFC(f),['SF=' sfNames{f} 'Ori=90' ]);
            end
        
            xlabel(plotHandlesSFC(f),'Frequency (Hz)');
            xlabel(plotHandlesSTA(f),'Time (ms)');
            xlabel(plotHandlesPAC5(1),'Phase Frequency from LFP data (Hz)');
        
            axis(plotHandlesFFC(f),[0 ffcFreq(end) 0 1]);
            axis(plotHandlesSFC(f),[0 ffcFreq(end) 0 1]);
            axis(plotHandlesSTA(f),[xsSTA(1) xsSTA(end) -80 40]);
    
            clim(plotHandlesPAC1(f),[min(dataPlot(:,1,:,:),[],'all') max(dataPlot(:,1,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC1(6));
            cb.Ticks = [min(dataPlot(:,1,:,:),[],'all') max(dataPlot(:,1,:,:),[],'all')];
            ylabel(cb,'MVLMI','FontWeight','bold','FontSize',8,'Rotation',270)
    
            clim(plotHandlesPAC2(f),[min(dataPlot(:,2,:,:),[],'all') max(dataPlot(:,2,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC2(6));
            cb.Ticks = [min(dataPlot(:,2,:,:),[],'all') max(dataPlot(:,2,:,:),[],'all')];
            ylabel(cb,'PLV','FontWeight','bold','FontSize',8,'Rotation',270)
    
            clim(plotHandlesPAC3(f),[min(dataPlot(:,3,:,:),[],'all') max(dataPlot(:,3,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC3(6));
            cb.Ticks = [min(dataPlot(:,3,:,:),[],'all') max(dataPlot(:,3,:,:),[],'all')];
            ylabel(cb,'GLM','FontWeight','bold','FontSize',8,'Rotation',270)
    
    
            clim(plotHandlesPAC4(f),[min(dataPlot(:,4,:,:),[],'all') max(dataPlot(:,4,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC4(6));
            cb.Ticks = [min(dataPlot(:,4,:,:),[],'all') max(dataPlot(:,4,:,:),[],'all')];
            ylabel(cb,'KLMI','FontWeight','bold','FontSize',8,'Rotation',270)
    
    
            clim(plotHandlesPAC5(f),[min(dataPlot(:,5,:,:),[],'all') max(dataPlot(:,5,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC5(6));
            cb.Ticks = [min(dataPlot(:,5,:,:),[],'all') max(dataPlot(:,5,:,:),[],'all')];
            ylabel(cb,'n-MVL','FontWeight','bold','FontSize',8,'Rotation',270)
            end

    elseif strcmp(tuning,'Size')
            sigNames={'0.1', '0.2', '0.4', '0.8', '1.6', '3.2', 'all'};
            for nsig=1:numOrientations
              for i=1:numCombinations
                
                plot(plotHandlesFFC(nsig),ffcFreq,ffcData{i}(nsig,:),'color',distanceColorNames(i));
                hold(plotHandlesFFC(nsig),'on');
        
                plot(plotHandlesSFC(nsig),sfcFreq,sfcData{i}(nsig,:),'color',distanceColorNames(i));
                hold(plotHandlesSFC(nsig),'on');
        
                plot(plotHandlesSTA(nsig),xsSTA,staData{i}(nsig,:),'color',distanceColorNames(i));
                hold(plotHandlesSTA(nsig),'on');
                
                contourf(plotHandlesPAC1(nsig),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(nsig,1,:,:)),30,'lines','none');  
                hold(plotHandlesPAC1(nsig),'on');

                contourf(plotHandlesPAC2(nsig),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(nsig,2,:,:)),30,'lines','none');  
                hold(plotHandlesPAC2(nsig),'on');

                contourf(plotHandlesPAC3(nsig),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(nsig,3,:,:)),30,'lines','none');  
                hold(plotHandlesPAC3(nsig),'on');

                contourf(plotHandlesPAC4(nsig),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(nsig,4,:,:)),30,'lines','none');  
                hold(plotHandlesPAC4(nsig),'on');

                contourf(plotHandlesPAC5(nsig),PhaseFreqVector, AmpFreqVector, squeeze(dataPlot(nsig,5,:,:)),30,'lines','none');  
                hold(plotHandlesPAC5(nsig),'on');
    
            end
        
            set(plotHandlesFFC(nsig),'XTickLabel',[]);
            if nsig==1
                ylabel(plotHandlesFFC(nsig),'FFC');
                ylabel(plotHandlesSFC(nsig),'SFC');
                ylabel(plotHandlesSTA(nsig),'STA');           
                
                ylabel(plotHandlesPAC5(nsig),sprintf('Amplitude Frequency (Hz)'));
                title(plotHandlesFFC(nsig),['Sigma=' sigNames{nsig} 'Ori=90' ]);
    
            else
                set(plotHandlesFFC(nsig),'YTickLabel',[]);
                set(plotHandlesSFC(nsig),'YTickLabel',[]);
                set(plotHandlesSTA(nsig),'YTickLabel',[]);
                set(plotHandlesPAC1(nsig),'YTickLabel',[]);
                set(plotHandlesPAC2(nsig),'YTickLabel',[]);
                set(plotHandlesPAC3(nsig),'YTickLabel',[]);
                set(plotHandlesPAC4(nsig),'YTickLabel',[]);
                set(plotHandlesPAC5(nsig),'YTickLabel',[]);
                title(plotHandlesFFC(nsig),['Sigma=' sigNames{nsig} 'Ori=90' ]);
            end
        
            xlabel(plotHandlesSFC(nsig),'Frequency (Hz)');
            xlabel(plotHandlesSTA(nsig),'Time (ms)');
            xlabel(plotHandlesPAC5(1),'Phase Frequency from LFP data (Hz)');
        
            axis(plotHandlesFFC(nsig),[0 ffcFreq(end) 0 1]);
            axis(plotHandlesSFC(nsig),[0 ffcFreq(end) 0 1]);
            axis(plotHandlesSTA(nsig),[xsSTA(1) xsSTA(end) -80 40]);
    
            clim(plotHandlesPAC1(nsig),[min(dataPlot(:,1,:,:),[],'all') max(dataPlot(:,1,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC1(7));
            cb.Ticks = [min(dataPlot(:,1,:,:),[],'all') max(dataPlot(:,1,:,:),[],'all')];
            ylabel(cb,'MVLMI','FontWeight','bold','FontSize',8,'Rotation',270)
    
            clim(plotHandlesPAC2(nsig),[min(dataPlot(:,2,:,:),[],'all') max(dataPlot(:,2,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC2(7));
            cb.Ticks = [min(dataPlot(:,2,:,:),[],'all') max(dataPlot(:,2,:,:),[],'all')];
            ylabel(cb,'PLV','FontWeight','bold','FontSize',8,'Rotation',270)
    
            clim(plotHandlesPAC3(nsig),[min(dataPlot(:,3,:,:),[],'all') max(dataPlot(:,3,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC3(7));
            cb.Ticks = [min(dataPlot(:,3,:,:),[],'all') max(dataPlot(:,3,:,:),[],'all')];
            ylabel(cb,'GLM','FontWeight','bold','FontSize',8,'Rotation',270)
    
    
            clim(plotHandlesPAC4(nsig),[min(dataPlot(:,4,:,:),[],'all') max(dataPlot(:,4,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC4(7));
            cb.Ticks = [min(dataPlot(:,4,:,:),[],'all') max(dataPlot(:,4,:,:),[],'all')];
            ylabel(cb,'KLMI','FontWeight','bold','FontSize',8,'Rotation',270)
    
    
            clim(plotHandlesPAC5(nsig),[min(dataPlot(:,5,:,:),[],'all') max(dataPlot(:,5,:,:),[],'all')]);
            cb=colorbar(plotHandlesPAC5(7));
            cb.Ticks = [min(dataPlot(:,5,:,:),[],'all') max(dataPlot(:,5,:,:),[],'all')];
            ylabel(cb,'n-MVL','FontWeight','bold','FontSize',8,'Rotation',270)
            end

           end
        end
end

end










