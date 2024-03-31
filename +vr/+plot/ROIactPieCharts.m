 %% 

% clear
function f = ROIactPieCharts(sData,sDataDir) 
% Select and load sData file

if nargin == 2
    sessionID = sData.sessionInfo.sessionID;
else
    [sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
    load(fullfile(sDataDir,sessionID));
    sessionID = sData.sessionInfo.sessionID;
end



%% 
%{
samplingRate = sData.daqdata.meta.samplingRate;
binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
homeBoxLength = abs(sData.stats.sessionAvs(1).plotXAxis(1) - binSize);
rewardZone = sData.behavior.rewardZone;
corridorLength = binNumber*binSize;
viewDistance = sData.behavior.viewDistance;

nHomeBoxBins = homeBoxLength/binSize;
binnedPosition = discretize(sData.daqdata.unityPosition,-homeBoxLength:binSize:-homeBoxLength+corridorLength);
trialStartIndexes = find(diff(sData.daqdata.unityPosition)< -10)+1;
allTrials = numel(trialStartIndexes)-1;
velocity = sData.behavior.signals.velocity;
%}

for fov = 1:1:length(sData.imdata)

    if isfield(sData.imdata(fov),'fovLocation')
        fovName = sData.imdata(fov).fovLocation;
    else
        fovName = '';
    end
    
    if length(sData.imdata) > 1
        fovName = ['Fov' num2str(fov) '-' fovName];
    end

    
nROIs = sData.imdata.nROIs;
nActiveROIs = numel(sData.imdata.activeROIs);
activeROIs = sData.imdata.activeROIs;

nInactiveROIs = numel(sData.imdata.inactiveROIs);
inactiveROIs = sData.imdata.inactiveROIs;

nTunedROIs = numel(sData.imdata.tunedROIs);
tunedROIs = sData.imdata.tunedROIs;

nUntunedROIs = numel(sData.imdata.untunedROIs);
untunedROIs = sData.imdata.untunedROIs;

%{
data = NaN;
data(1) = numel(intersect(inactiveROIs,untunedROIs));
data(2) = numel(intersect(inactiveROIs,tunedROIs));
data(3) = numel(intersect(activeROIs,tunedROIs));
data(4) = numel(intersect(activeROIs,untunedROIs));

numel(intersect(find(activeTrialFractionCtrl > threshold),find(peakFreqsCtrl > freqThreshold)))

label = {'Inactive untuned','Inactive tuned','Active tuned','Active untuned'};

figure
subplot(2,2,1)
pie(data,[1 0 0 0],label)
subplot(2,2,3)
pie(data,[0 1 0 0],label)
subplot(2,2,4)
pie(data,[0 0 1 0],label)
subplot(2,2,2)
pie(data,[0 0 0 1],label)
%}


f = figure('Color','white','Position',[0 0 600 400]);
subplot(1,2,1)
data = NaN;
data(1) = nInactiveROIs;
data(2) = nActiveROIs;
label = {'Inactive','Active'};
pie(data,label)
% title(['All ROIs: ' num2str(nROIs) newline 'Active: ' num2str(nActiveROIs) newline 'Inactive: ' num2str(nInactiveROIs)])
title(['Active ROIs: ' num2str(nActiveROIs) ' / ' num2str(nROIs) newline 'Inactive ROIs: ' num2str(nInactiveROIs) ' / ' num2str(nROIs)])


subplot(1,2,2)
data = NaN;
data(1) = numel(intersect(inactiveROIs,untunedROIs));
data(2) = numel(intersect(inactiveROIs,tunedROIs));
data(3) = numel(intersect(activeROIs,tunedROIs));
data(4) = numel(intersect(activeROIs,untunedROIs));
label = {'Inactive untuned','Inactive tuned','Active tuned','Active untuned'};
%label = {['Inactive' newline 'untuned'],['Inactive' newline 'tuned'],['Active' newline 'tuned'],['Active' newline 'untuned']};
pie(data,label)
% title(['All ROIs: ' num2str(nROIs) newline 'Tuned: ' num2str(nTunedROIs) newline 'Untuned: ' num2str(nUntunedROIs)])
title(['Tuned ROIs: ' num2str(nTunedROIs) ' / ' num2str(nROIs) newline 'Untuned ROIs: ' num2str(nUntunedROIs) ' / ' num2str(nROIs)])

suptitle([sessionID(1:17) ' - ' fovName])



saveas(f,strcat(fullfile(sDataDir,[sessionID(1:17) ' - ' fovName '-ROIactPieCharts']), '.png'));
%close(gcf)
    
end

end

