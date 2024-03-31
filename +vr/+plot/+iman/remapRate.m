function sData = remapRate(sData,sDataDir)

for fov = 1:1:length(sData.imdata)

    blockData = sData.imdata(fov).placeCells;
    %sessionID = sData.sessionInfo.sessionID;
    
    
    
%% Fig 1 RemapRate 2 Blocks
figure('Color','white','Position',[0 0 400 300])

hold on
for t = 1:2
    indTrials = blockData(t).inductionTrial;
    %indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos <= 15)' blockData(t).inductionTrial(blockData(t).peakPos >= 90)']'; % identical
    %indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos > 15)' blockData(t).inductionTrial(blockData(t).peakPos < 90)']'; % different    
    indTrials = indTrials(indTrials > 0); % exclude zeros
    
    x = 0:1:max(indTrials);
    y = zeros(size(x));
    for i = 1:numel(x)
        y(i) = sum(indTrials <= x(i));
    end
    y = y/numel(indTrials);
    stairs(x,y,'-','LineWidth',1.5)
end
xlim([0 30])
xlabel('Place field onset trial')
ylabel('Cummulative fraction')
title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
legend({sData.trials.contextsMeta.name})

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation '_mappingRate2Blocks']),'.png'));

%% Fig 2 RemapRate 4 Blocks
figure('Color','white','Position',[0 0 400 300])

hold on
for t = 1:4
    indTrials = blockData(t).inductionTrial;
    %indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos <= 15)' blockData(t).inductionTrial(blockData(t).peakPos >= 90)']'; % identical
    %indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos > 15)' blockData(t).inductionTrial(blockData(t).peakPos < 90)']'; % different    
    indTrials = indTrials(indTrials > 0); % exclude zeros
    
    x = 0:1:max(indTrials);
    y = zeros(size(x));
    for i = 1:numel(x)
        y(i) = sum(indTrials <= x(i));
    end
    y = y/numel(indTrials);
    stairs(x,y,'-','LineWidth',1.5)
end
xlim([0 30])
xlabel('Place field onset trial')
ylabel('Cummulative fraction')
title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
legend('Context A1','Context B1','Context A2','Context B2')

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation '_mappingRate4Blocks']),'.png'));

%% Fig 3 RemapRate 2 Blocks - Visually identical part of the track
figure('Color','white','Position',[0 0 400 300])

hold on
for t = 1:2
    %indTrials = blockData(t).inductionTrial;
    indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos <= 15)' blockData(t).inductionTrial(blockData(t).peakPos >= 90)']'; % identical
    %indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos > 15)' blockData(t).inductionTrial(blockData(t).peakPos < 90)']'; % different    
    indTrials = indTrials(indTrials > 0); % exclude zeros
    
    x = 0:1:max(indTrials);
    y = zeros(size(x));
    for i = 1:numel(x)
        y(i) = sum(indTrials <= x(i));
    end
    y = y/numel(indTrials);
    stairs(x,y,'-','LineWidth',1.5)
end
xlim([0 30])
xlabel('Place field onset trial')
ylabel('Cummulative fraction')
title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
legend('Context A1','Context B1','Context A2','Context B2')

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation '_mappingRateIdenticalPart']),'.png'));

%% Fig 4 RemapRate 2 Blocks - Visually different part of the track
figure('Color','white','Position',[0 0 400 300])

hold on
for t = 1:2
    %indTrials = blockData(t).inductionTrial;
    %indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos <= 15)' blockData(t).inductionTrial(blockData(t).peakPos >= 90)']'; % identical
    indTrials =  [blockData(t).inductionTrial(blockData(t).peakPos > 15)' blockData(t).inductionTrial(blockData(t).peakPos < 90)']'; % different    
    indTrials = indTrials(indTrials > 0); % exclude zeros
    
    x = 0:1:max(indTrials);
    y = zeros(size(x));
    for i = 1:numel(x)
        y(i) = sum(indTrials <= x(i));
    end
    y = y/numel(indTrials);
    stairs(x,y,'-','LineWidth',1.5)
end
xlim([0 30])
xlabel('Place field onset trial')
ylabel('Cummulative fraction')
title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
legend('Context A1','Context B1','Context A2','Context B2')

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation '_mappingRateDifferentPart']),'.png'));



%%
% extract data for session comparision plot
%{
Diffremap(4).name = sessionID(1:5);
Diffremap(4).x = x;
Diffremap(4).y = y;



Idremap(1).name = sessionID(1:5);
Idremap(1).x = x;
Idremap(1).y = y;

Idremap(2).name = sessionID(1:5);
Idremap(2).x = x;
Idremap(2).y = y;

Aremap(4).name = sessionID(1:5);
Aremap(4).x = x;
Aremap(4).y = y;
%}

%% Plot place field numbers and plase cell fractions
nROIs = sData.imdata(fov).nROIs;

pieData = NaN;
for trialType = 1:1:size(sData.imdata(fov).placeCells,2)
    placeCells = [];
    nFields = [];
    for roi = 1:1:sData.imdata(fov).nROIs
        if sData.imdata(fov).roiMeta(roi).placeCell(trialType)
            placeCells = [placeCells roi];
            nFields = [nFields sData.imdata(fov).roiMeta(roi).nFields(trialType)];
        end
    end
    blockData(trialType).placeCells = placeCells;
    blockData(trialType).nFields = nFields;
    blockData(trialType).nPlaceCells = numel(placeCells);
    blockData(trialType).placeCellFraction = blockData(trialType).nPlaceCells/nROIs;
    
    pieData(trialType,1) = nROIs - numel(placeCells); % 0 fields
    pieData(trialType,2) = numel(nFields(nFields == 1));
    pieData(trialType,3) = numel(nFields(nFields == 2));
    pieData(trialType,4) = numel(nFields(nFields > 2));
    
    blockData(trialType).fieldCounts012more = pieData(trialType,:);
end





figure('Color','white','Position',[0 0 400 300])
x = categorical({'1 - A','2 - B','3 - A','4 - B'});
bar(x,pieData(:,2:4)/nROIs,'stacked')
ylim([0 1])
ylabel('Fraction of ROIs')
xlabel('Trial blocks')
legend({'1 Place fields','2 Place fields','>2 Place fields'},'Location','northeast')
title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation '_placeCellCountsBar']),'.png'));




sData.imdata(fov).placeCells = blockData;




end


end