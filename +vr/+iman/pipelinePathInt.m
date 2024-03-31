%% Image analysis pipeline for Path integration VR task

 %% This version was used for the "random homebox" task (2020 winter)

clear

% Select and load sData file
[sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
load(fullfile(sDataDir,sessionID));
sessionID = strsplit(sessionID,'.mat');
sessionID = sessionID{1};

% vr.plot.iman.meanOptoEffectOnPCs(sData,sDataDir)

% sData = vr.plot.iman.remapPlotsRoiMatricesOpto(sData,sDataDir);
% save(fullfile(sDataDir,sessionID),'sData');  close all




sData = vr.downsampleBehavior(sData); % save(fullfile(sDataDir,sessionID),'sData');

%% quality check

sData = vr.iman.imQualityChecks(sData, sDataDir);


% save(fullfile(sDataDir,sessionID),'sData');
sData = vr.iman.createRoiMatricesContexts(sData);
% save(fullfile(sDataDir,sessionID),'sData');

%% Correction
% sData = vr.iman.createRoiMatricesPathInt(sData);
%{
if isfield(sData.trials.contextsMeta,'placeCells')
    sData.trials.contextsMeta = rmfield(sData.trials.contextsMeta,'placeCells');
end
if isfield(sData.trials.contextsMeta,'nPlaceCells')
    sData.trials.contextsMeta = rmfield(sData.trials.contextsMeta,'nPlaceCells');
end
if isfield(sData.trials.contextsMeta,'placeCellFraction')
    sData.trials.contextsMeta = rmfield(sData.trials.contextsMeta,'placeCellFraction');
end
%}


%% Add FOV location
% fovLocation = {'PPC', 'M2(RSC)', 'RSC', 'RSC'};
% fovLocation = {'PPC', 'V2', 'RSC', 'M2'};
% fovLocation = {'RSC', 'V2', 'M2(RSC)',  'PPC(V2)'};
% fovLocation = {'V2(PPC)', 'PPC', 'RSC', 'M2'};
% fovLocation = {'V2', 'RSC', 'M2',  'PPC-V2'};
% fovLocation = {'V2(RSC)', 'PPC', 'RSC', 'M2'};
% fovLocation = {'V2', 'PPC', 'RSC', 'M2'};
% fovLocation = {'HPC'};
% fovLocation = {'RSC'};
% fovLocation = {'V1', 'V2', 'RSC', 'PPC'};
% fovLocation = {'V2', 'V2', 'RSC', 'M2'};
% fovLocation = {'V2', 'V2', 'PPC', 'RSC'};
% fovLocation = {'RSC(antL)', 'RSC(antR)', 'RSC(postL)', 'RSC(postR)'};
% fovLocation = {'RSC(antL)', 'RSC(postL)', 'RSC(antR)', 'RSC(postR)'};

for fov = 1:1:length(sData.imdata)
    sData.imdata(fov).fovLocation = fovLocation{fov};
end


% save(fullfile(sDataDir,sessionID),'sData');



sData = vr.iman.classifyROIs(sData, sDataDir);
%save(fullfile(sDataDir,sessionID),'sData');
close all


sData = vr.plot.iman.remapRate(sData,sDataDir);
% sData = vr.plot.iman.remapRateOpto(sData,sDataDir);
%save(fullfile(sDataDir,sessionID),'sData');

vr.plot.iman.individualROIs4TrialBlocks(sData,sDataDir)
% vr.plot.iman.individualROIs5TrialBlocksOpto(sData,sDataDir); vr.plot.iman.meanOptoEffectOnPCs(sData,sDataDir)

%vr.plot.iman.trialToTrialSimilarity
sData = vr.plot.iman.trialToTrialSimilarity(sData,sDataDir);
save(fullfile(sDataDir,sessionID),'sData');




% vr.plot.ROIactPieCharts(sData,sDataDir) % old version with tuning
vr.plot.iman.ROIactPieChart(sData,sDataDir);

if sData.behavior.trialMatrices.meta.binNumber > 125
vr.plot.iman.remapPlotsRoiMatricesLong(sData,sDataDir) 
else
vr.plot.iman.remapPlotsRoiMatrices(sData,sDataDir) 
end

% sData = vr.plot.iman.remapPlotsRoiMatricesOpto(sData,sDataDir);
save(fullfile(sDataDir,sessionID),'sData');

close all





figure; for i = 1:20;  subplot(4,5,i); imagesc(sData.imdata.binnedRoisDff(:,:,i)); end






% plot DFF signals
nFOVs = length(sData.imdata);

for f = 1:1:nFOVs
    [~, ROIsSortedSNR] = sort(sData.imdata(f).roiStat.signalToNoise);
    %sDataSNRSort.imdata.roiSignals(2).dff; sData.imdata.roiSignals(2).dff(ROIsSortedSNR,:);
    range = 1:4500; %5001:10000;
    
    for i = 0:1:floor(sData.imdata(f).nROIs/10)
        if i*10+10 < sData.imdata(f).nROIs
            subset = i*10+1:1:i*10+10;
        else
            subset = i*10+1:1:sData.imdata(f).nROIs;
        end
        sDataTemp.imdata = sData.imdata(f);
        sDataTemp.sessionInfo = sData.sessionInfo;
        try
            vr.plot.dffSignals(sDataTemp,ROIsSortedSNR(subset),range);
        catch
        end
        if ~isfolder(fullfile(sDataDir,['dffSignals_' num2str(f)]))
            mkdir(fullfile(sDataDir,['dffSignals_' num2str(f)]));
        end
        saveas(gcf,strcat(fullfile(sDataDir,['dffSignals_' num2str(f)],['dffSignals-' num2str(i)]),'.png'));
        close gcf
    end
    
end










































