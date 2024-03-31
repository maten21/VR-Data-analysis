function f = individualROIs4TrialBlocks(sData,sDataDir) 



for fov = 1:1:length(sData.imdata)
%% Plot individual ROIs 4 trial blocks Remapping
% subsetOfInterest = 1:1:nROIs;

% Place cell in either block:
% subsetOfInterest = union(union(sData.imdata(fov).placeCells(1).placeCells,sData.imdata(fov).placeCells(2).placeCells),...
%    union(sData.imdata(fov).placeCells(3).placeCells,sData.imdata(fov).placeCells(4).placeCells));


try
            %A =      [sData.imdata(fov).roiMeta.identPartCorrCoefA];
            %B =      [sData.imdata(fov).roiMeta.identPartCorrCoefB];
            %AB =     [sData.imdata(fov).roiMeta.identPartCorrCoefAB];
            ASign =  [sData.imdata(fov).roiMeta.identPartIsSignCorrA];
            BSign =  [sData.imdata(fov).roiMeta.identPartIsSignCorrB];
            ABSign = [sData.imdata(fov).roiMeta.identPartIsSignCorrAB];
        %end
        
        
            ASignRois = find(ASign);
            BSignRois = find(BSign);
            ABSignRois = find(ABSign);
            nROIs = sData.imdata(fov).nROIs;
            
            remapRoisA = setdiff(ASignRois,ABSignRois); 
            %nonRemapRoisA = intersect(ASignRois,ABSignRois);
            remapRoisB = setdiff(BSignRois,ABSignRois); 
            %nonRemapRoisB = intersect(BSignRois,ABSignRois);

            remapRois = union(remapRoisA,remapRoisB);
catch
end





 subsetOfInterest = sData.imdata(fov).activeROIs;
 tag = 'activeCells_dff';

%subsetOfInterest = remapRois;
%tag = 'remapRoisIdent_dff';


smoothSpan = 5;

binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;

rew = nanmean(sData.behavior.trialMatrices.rewardInBin);
rewardZones = sData.stats.sessionAvs(1).plotXAxis(find(rew));
rewardZones = rewardZones - binSize;

rewardZones(rewardZones==0) = binNumber*binSize;



if isfield(sData.imdata(fov),'fovLocation')
    fovLoc = ['_' sData.imdata(fov).fovLocation];
else
    fovLoc = '';
end

for roi = 1:1:numel(subsetOfInterest)
    %data = permute(binnedRoisLickAlignedDff,[3 1 2]);
    figure('Color','white','Position',[0 0 800 800]);
      
    
    data = smoothdata(sData.imdata(fov).binnedRoisDff(:,:,subsetOfInterest(roi)),2,'gaussian',smoothSpan);
    clims = quantile(data(:),[0.02 0.98]);
    scaling = clims(2)-clims(1);
    Xax = sData.stats.sessionAvs(1).plotXAxis;
    
    subplot(2,2,1)
    
    trials = sData.trials.contextsMeta(1).trials;
    nTrials = numel(trials);    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan)/scaling *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)

    if length(rewardZones) > 2
        rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(3),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    else
        rectangle('Position',[rewardZones(1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    end
    
    %xlim([-320 180])
    ylim([1 nTrials])
    
    caxis(clims)
    %text(140,-5,['\fontsize{13}' num2str(clims(2),2)])
    % text(-120,-5,['\fontsize{14}' 'Control'])
   % xlabel('\fontsize{13}Track position (cm)');
    ylabel('\fontsize{13}Trials');
    t = title('Block 1 - Context A');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    subplot(2,2,2)
    
    trials = sData.trials.contextsMeta(2).trials;
    nTrials = numel(trials);    
    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan)/scaling *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)

    if length(rewardZones) > 2
        rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        %rectangle('Position',[rewardZones(4),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(3),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    else
        rectangle('Position',[rewardZones(1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    end
    
    %xlim([-320 180])
    ylim([1 nTrials])

    caxis(clims)
    text(250,0,['\fontsize{13}' num2str(clims(2),2)])
    %text(-120,-5,['\fontsize{14}' 'Control'])
   % xlabel('\fontsize{13}Track position (cm)');
   % ylabel('\fontsize{13}Trials');
    t = title('Block 2 - Context B');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    
    subplot(2,2,3)
    
    trials = sData.trials.contextsMeta(3).trials;
    nTrials = numel(trials);    
    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan)/scaling *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)

    if length(rewardZones) > 2
        rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(3),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    else
        rectangle('Position',[rewardZones(1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    end
    
    %xlim([-320 180])
    ylim([1 nTrials])
    
    caxis(clims)
    % text(140,-5,['\fontsize{13}' num2str(clims(2),2)])
    %text(-120,-5,['\fontsize{14}' 'Control'])
    xlabel('\fontsize{13}Track position (cm)');
    ylabel('\fontsize{13}Trials');
    t = title('Block 3 - Context A');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    
    subplot(2,2,4)
    
    trials = sData.trials.contextsMeta(4).trials;
    nTrials = numel(trials);    
    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan)/scaling *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)

    if length(rewardZones) > 2
        rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        %rectangle('Position',[rewardZones(4),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(3),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    else
        rectangle('Position',[rewardZones(1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
        rectangle('Position',[rewardZones(2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    end
    
    %xlim([-320 180])
    ylim([1 nTrials])

    caxis(clims)
    % text(140,-5,['\fontsize{13}' num2str(clims(2),2)])
    %text(-120,-5,['\fontsize{14}' 'Control'])
    xlabel('\fontsize{13}Track position (cm)');
   % ylabel('\fontsize{13}Trials');
    t = title('Block 4 - Context B');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    suptitle(['ROI: ' num2str(subsetOfInterest(roi))])
 
    
    
    

    if ~exist([sDataDir 'ROIsIn4Blocks_' 'Fov' num2str(fov) fovLoc '_' tag], 'dir')
        mkdir([sDataDir 'ROIsIn4Blocks_' 'Fov' num2str(fov) fovLoc '_' tag]);
    end    
    saveas(gcf,strcat(fullfile(sDataDir,['ROIsIn4Blocks_' 'Fov' num2str(fov) fovLoc '_' tag],['ROI-' num2str(subsetOfInterest(roi))]), '.png'));
    close(gcf)
    
    
    
    
end







end





