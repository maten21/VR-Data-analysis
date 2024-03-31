
% Structure containing the metadata related to the contexts trial blocks: 
% sData.trials.contextsMeta
%
% Structure containing binned behavior data: 
% sData.behavior.trialMatrices
%    sData.behavior.trialMatrices.binVel           % rinning speed
%    sData.behavior.trialMatrices.lickFreqInBin    % lick rate
%
% Structure containing all imaging related data: 
% sData.imdata
%    Data.imdata.binnedRoisDff                     % binned 3D ROI matrix, with the dimensions: sData.imdata.binnedRoisDff(trials,trackPosBin,rois)



%% Example plot 1:
% Plot the first ROI activity in the first trial block

blockIndex = 1; % nTrialBlocks = length(sData.trials.contextsMeta);
trials = sData.trials.contextsMeta(blockIndex).trials;
trackPosBin = 1:1:sData.imdata.binNumber;
roiIndex = 1; %

roiMatrix = sData.imdata.binnedRoisDff(trials,trackPosBin,roiIndex);

% For the plot X axis generate an array with the true position in cm instead of the bin number
binSize = sData.imdata.binSize;
nBins = sData.imdata.binNumber;

trackPosCm = binSize:binSize:nBins*binSize;


figure
imagesc(trackPosCm,trials,roiMatrix)


%% Example plot 2:
% Plot the first ROI activity in all trial blocks with subplots

nTrialBlocks = length(sData.trials.contextsMeta); % number of thial blocks in the experiment

% For the plot X axis generate an array with the true position in cm instead of the bin number
binSize = sData.imdata.binSize;
nBins = sData.imdata.binNumber;
trackPosCm = binSize:binSize:nBins*binSize;
roiIndex = 1; %

subplotRows = 2;
subplotCol = 2;

figure
for i = 1:1:nTrialBlocks
subplot(subplotRows,subplotCol,i)

trials = sData.trials.contextsMeta(i).trials;

roiMatrix = sData.imdata.binnedRoisDff(trials,:,roiIndex);

imagesc(trackPosCm,trials,roiMatrix)

title(sData.trials.contextsMeta(i).name)
end


