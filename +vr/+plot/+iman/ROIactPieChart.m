% Plot active inactive rois pie charts for each FOV
% clear
function f = ROIactPieChart(sData,sDataDir) 
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
nFovs = length(sData.imdata);

for fov = 1:1:nFovs

    if isfield(sData.imdata(fov),'fovLocation')
        fovName = sData.imdata(fov).fovLocation;
    else
        fovName = '';
    end
    
    if length(sData.imdata) > 1
        fovName = ['Fov' num2str(fov) '-' fovName];
    end

    
nROIs = sData.imdata(fov).nROIs;
nActiveROIs = numel(sData.imdata(fov).activeROIs);
activeROIs = sData.imdata(fov).activeROIs;

nInactiveROIs = numel(sData.imdata(fov).inactiveROIs);
inactiveROIs = sData.imdata(fov).inactiveROIs;

%nTunedROIs = numel(sData.imdata(fov).tunedROIs);
%tunedROIs = sData.imdata(fov).tunedROIs;

%nUntunedROIs = numel(sData.imdata(fov).untunedROIs);
%untunedROIs = sData.imdata(fov).untunedROIs;

%data = NaN;
data(fov,1) = nInactiveROIs;
data(fov,2) = nActiveROIs;
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
end


f = figure('Color','white','Position',[0 0 300*nFovs 400]);

for fov = 1:1:nFovs
subplot(1,nFovs,fov)

label = {'Inactive','Active'};
pie(data(fov,:),label)
% title(['All ROIs: ' num2str(nROIs) newline 'Active: ' num2str(nActiveROIs) newline 'Inactive: ' num2str(nInactiveROIs)])
title(['Fov' num2str(fov) ' - ' sData.imdata(fov).fovLocation newline ... 
    'Active ROIs: ' num2str(data(fov,2)) ' / ' num2str(num2str(sData.imdata(fov).nROIs)) newline ... 
    'Inactive ROIs: ' num2str(data(fov,1)) ' / ' num2str(sData.imdata(fov).nROIs)])

%{
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
%}

end
suptitle([sessionID(1:17) ])



saveas(f,strcat(fullfile(sDataDir,[sessionID(1:17)  '-ROIactInFovsPieCharts']), '.png'));
%close(gcf)
    


end
