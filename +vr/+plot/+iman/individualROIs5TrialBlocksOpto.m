function f = individualROIs5TrialBlocksOpto(sData,sDataDir) 



for fov = 1:1:length(sData.imdata)
%% Plot individual ROIs 5 trial blocks Remapping with opto
% subsetOfInterest = 1:1:nROIs;

% Place cell in either block:
% subsetOfInterest = union(union(sData.imdata(fov).placeCells(1).placeCells,sData.imdata(fov).placeCells(2).placeCells),...
%    union(sData.imdata(fov).placeCells(3).placeCells,sData.imdata(fov).placeCells(4).placeCells));

subsetOfInterest = sData.imdata(fov).activeROIs;

smoothSpan = 1;

binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
rewardZones = [];

for context = 1:1:length(sData.trials.contextsMeta)
rew = mean(sData.behavior.trialMatrices.rewardInBin(sData.trials.contextsMeta(context).trials,:));
rew(isnan(rew)) = 0;
rewardZ = sData.stats.sessionAvs(1).plotXAxis(find(rew));

rewardZones(context,1) = rewardZ(1); 
rewardZones(context,2) = rewardZ(find(diff(discretize(rewardZ,2)))+1); 


end
rewardZones = rewardZones - binSize;


tag = 'activeCells_dff';
if isfield(sData.imdata(fov),'fovLocation')
    fovLoc = ['_' sData.imdata(fov).fovLocation];
else
    fovLoc = '';
end

for roi = 1:1:numel(subsetOfInterest)
    %data = permute(binnedRoisLickAlignedDff,[3 1 2]);
    figure('Color','white','Position',[0 0 1200 800]);
      
    
    data = sData.imdata(fov).binnedRoisDff(:,:,subsetOfInterest(roi));
    clims = quantile(data(:),[0 0.99]);
    Xax = sData.stats.sessionAvs(1).plotXAxis;
    
    subplot(2,3,2)
    
    trials = sData.trials.contextsMeta(1).trials;
    nTrials = numel(trials);    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan) *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)
    % rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    % rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(1,1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(1,2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');

    %xlim([-320 180])
    ylim([1 nTrials])
    
    caxis(clims)
    %text(140,-5,['\fontsize{13}' num2str(clims(2),2)])
    % text(-120,-5,['\fontsize{14}' 'Control'])
   % xlabel('\fontsize{13}Track position (cm)');
    ylabel('\fontsize{13}Trials');
    t = title('Block 1 - Fam - No-stim.');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    
    subplot(2,3,3)
    
    trials = sData.trials.contextsMeta(5).trials;
    nTrials = numel(trials);    
    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan) *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)
    % rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    % rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(5,1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    %rectangle('Position',[rewardZones(4),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(5,2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    
    %xlim([-320 180])
    ylim([1 nTrials])

    caxis(clims)
    text(250,0,['\fontsize{13}' num2str(clims(2),2)])
    %text(-120,-5,['\fontsize{14}' 'Control'])
   % xlabel('\fontsize{13}Track position (cm)');
   % ylabel('\fontsize{13}Trials');
    t = title('Block 5 - Fam - Opto','Color','r');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    
    subplot(2,3,4)
    
    trials = sData.trials.contextsMeta(2).trials;
    nTrials = numel(trials);    
    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan) *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)
    % rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    % rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(2,1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(2,2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');

    %xlim([-320 180])
    ylim([1 nTrials])
    
    caxis(clims)
    % text(140,-5,['\fontsize{13}' num2str(clims(2),2)])
    %text(-120,-5,['\fontsize{14}' 'Control'])
    xlabel('\fontsize{13}Track position (cm)');
    ylabel('\fontsize{13}Trials');
    t = title('Block 2 - New - Opto','Color','r');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    
    subplot(2,3,5)
    
    trials = sData.trials.contextsMeta(3).trials;
    nTrials = numel(trials);    
    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan) *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)
    % rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    % rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(3,1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(3,2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');

    %xlim([-320 180])
    ylim([1 nTrials])
    
    caxis(clims)
    % text(140,-5,['\fontsize{13}' num2str(clims(2),2)])
    %text(-120,-5,['\fontsize{14}' 'Control'])
    xlabel('\fontsize{13}Track position (cm)');
    ylabel('\fontsize{13}Trials');
    t = title('Block 3 - New - No-stim');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    
    subplot(2,3,6)
    
    trials = sData.trials.contextsMeta(4).trials;
    nTrials = numel(trials);    
    
    %Ydata = binnedRoisDeconvRate(:,:,sortedROIs(roi))/min(ciaDeconv(sortedROIs(roi),ciaDeconv(sortedROIs(roi),:)>0));
    Ydata = data(trials,:);
    imagesc(Xax,nTrials:1,Ydata)
    hold on
    plot(Xax,smoothdata(nanmean(Ydata),2,'gaussian',smoothSpan) *-nTrials/3+nTrials,'w-','LineWidth',1.5)
    %plot(Xax,okada(okada(nanmean(Ydata),2),2)*-nTrials/10+nTrials,'w-','LineWidth',1.5)
    % rectangle('Position',[50,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    % rectangle('Position',[190,0,2,nTrials],'FaceColor',[1 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(4,1),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    %rectangle('Position',[rewardZones(4),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    rectangle('Position',[rewardZones(4,2),0,2,nTrials],'FaceColor',[0 1 1],'EdgeColor','none');
    
    %xlim([-320 180])
    ylim([1 nTrials])

    caxis(clims)
    % text(140,-5,['\fontsize{13}' num2str(clims(2),2)])
    %text(-120,-5,['\fontsize{14}' 'Control'])
    xlabel('\fontsize{13}Track position (cm)');
   % ylabel('\fontsize{13}Trials');
    t = title('Block 4 - New - Opto','Color','r');
    t.FontWeight = 'normal';
    t.FontSize = 14;
    
    suptitle(['ROI: ' num2str(subsetOfInterest(roi))])
 
    
    
    

    if ~isfolder(fullfile(sDataDir, ['ROIsIn4Blocks_' 'Fov' num2str(fov) fovLoc '_' tag]))
        mkdir(fullfile(sDataDir, ['ROIsIn4Blocks_' 'Fov' num2str(fov) fovLoc '_' tag]));
    end    
    saveas(gcf,strcat(fullfile(sDataDir,['ROIsIn4Blocks_' 'Fov' num2str(fov) fovLoc '_' tag],['ROI-' num2str(subsetOfInterest(roi))]), '.png'));
    close(gcf)
    
    
    
    
end







end





