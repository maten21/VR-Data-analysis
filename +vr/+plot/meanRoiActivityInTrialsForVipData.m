

clear

% Select and load sData file
[sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
load(fullfile(sDataDir,sessionID));
sessionID = strsplit(sessionID,'.mat');
sessionID = sessionID{1};




blockStarts = [sData.trials.contextsMeta.blockStart];
nBlocTrials = [sData.trials.contextsMeta.nTrials];
Xlims = [0 sum([sData.trials.contextsMeta.nTrials])+1];

Xax = 1:1:sum([sData.trials.contextsMeta.nTrials]);

    
fig = figure('Color','white','Position',[0 0 400 300]);

fov = 1;
hold on
A = nanmean(permute(sData.imdata(fov).binnedRoisDff,[1,3,2]),3); % get rid of the spatial bins
rectangle('Position',[blockStarts(2),min(nanmean(A,2)),nBlocTrials(2),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
%rectangle('Position',[blockStarts(4),min(nanmean(A,2)),nBlocTrials(4),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
%plot(smoothdata(nanmean(A,2),'gaussian',5),'LineWidth',2)
plotMeanStd(Xax,A,lines(1),2)
xlim(Xlims)
ylim([min(nanmean(A,2))-0.01*min(nanmean(A,2)) max(nanmean(A,2))+0.01*min(nanmean(A,2))])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
%title(['FOV: ' num2str(fov) ' ' sData.imdata(fov).fovLocation])

    
title(sData.sessionInfo.sessionID(1:17))

saveas(fig,fullfile(sDataDir,['meanActivityAllROIs', '.png']));
    