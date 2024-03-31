% Create trial matrices for all ROIs. Part of image analysis pipeline.
% iFov is an optional input for multi FOV mesoscope recordings
% 2022.11.22
% 
function sData = createRoiMatricesPathInt(sData)

% SET VELOCITY THRESHOLD
velThreshold = 0.5; % Do not consider calcium signal if the mouse stops, or walks slower then this value.

%% Prepare data and create trial matrices for all ROIs

samplingRate = sData.daqdata.meta.samplingRate;
binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
homeBoxLength = abs(sData.stats.sessionAvs(1).plotXAxis(1) - binSize);
%rewardZone = sData.behavior.rewardZone;
corridorLength = binNumber*binSize;
%viewDistance = sData.behavior.viewDistance;

%nHomeBoxBins = homeBoxLength/binSize;
binnedPosition = discretize(sData.behavior.signals(1).corridorPosition,-homeBoxLength:binSize:-homeBoxLength+corridorLength);
trialStartIndexes = find(diff(sData.daqdata.unityPosition)< -10)+1;
allTrials = numel(trialStartIndexes)-1;
%velocity = sData.behavior.signals.velocity;
velocity = smoothdata(sData.behavior.signals(1).velocity,'gaussian',samplingRate/30); % Gaussian filter 30 Hz


%%

nFOVs =  length(sData.imdata);

for i = 1:1:nFOVs


% ROI signals 
dff = sData.imdata(i).roiSignals(2).dff;
ciaDeconv = sData.imdata(i).roiSignals(2).deconv;
deconvRate = sData.imdata(i).roiSignals(2).actRateDeconv;
nROIs = numel(dff(:,1));

   
try
    if isequal(sData.imdata(i).meta.imagingSystem(1:5),'Path1')
        scanFrameRate = sData.imdata(i).meta.scanFrameRate;
        samplePerFrame = sData.daqdata.meta.fs/scanFrameRate;
        frameIndexes = vr.fixFrameIndexes(sData.daqdata.frameIndex, samplePerFrame);
    elseif isequal(sData.imdata(i).meta.imagingSystem(1:5),'Path2')
        scanFrameRate = sData.imdata(i).meta.scanFrameRate;
        samplePerFrame = sData.daqdata.meta.fs/scanFrameRate;
        frameIndexes = vr.fixFrameIndexes(sData.daqdata.frameIndex2, samplePerFrame);
    end
catch
    scanFrameRate = sData.imdata.meta.fps;
    samplePerFrame = sData.daqdata.meta.fs/scanFrameRate;
    frameIndexes = vr.fixFrameIndexes(sData.daqdata.frameIndex, samplePerFrame);
end



frameSignalInd(1:numel(sData.daqdata.frameIndex)) = zeros; % This modified frame signal contains the index in the frame array ..001..0002.. 003

for ind = 1:1:numel(frameIndexes)

    frameSignalInd(frameIndexes(ind)) = ind;
    
end


binnedRoisDff = nan(allTrials,binNumber,nROIs);
binnedRoisDeconv = nan(allTrials,binNumber,nROIs);
binnedRoisDeconvRate = nan(allTrials,binNumber,nROIs);


% Create roi matrices DFF and Deconv
for r = 1:1:allTrials
    tempPos = binnedPosition(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    tempFrameSInd = frameSignalInd(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    tempVel = velocity(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    tempFrameSInd(find(tempVel < velThreshold)) = 0;
    tempBinPos = tempPos(find(tempFrameSInd > 1));
    
    for roi = 1:1:nROIs
        tempSignalDff = dff(roi,tempFrameSInd(tempFrameSInd>0));
        tempSignalDeconv = ciaDeconv(roi,tempFrameSInd(tempFrameSInd>0));
        tempSignalDeconvRate = deconvRate(roi,tempFrameSInd(tempFrameSInd>0));
        for c = 1:1:binNumber
            binnedRoisDff(r,c,roi) = mean(tempSignalDff(find(tempBinPos == c)));
            binnedRoisDeconv(r,c,roi) = mean(tempSignalDeconv(find(tempBinPos == c)));  
            binnedRoisDeconvRate(r,c,roi) = mean(tempSignalDeconvRate(find(tempBinPos == c)));  
        end
        clear('tempSignalDff','tempSignalDeconv','tempSignalDeconvRate');
    end
    clear('tempPos','tempFrameSInd','tempBinPos','tempVel');
end



% EXTEND TRIALS TO ALIGN DATA TO FIRST LICK
% 
fitstTrackBin = find(sData.stats.sessionAvs(1).plotXAxis == binSize);
lastTrackBin = find(sData.stats.sessionAvs(1).plotXAxis == sData.behavior.rewardZone);
firstLicks(1,allTrials) = NaN; 


for t = 1:1:allTrials
    firstLickPos = fitstTrackBin + min(find(sData.behavior.trialMatrices.licksInBin(t,fitstTrackBin:lastTrackBin)>0)) -1;
    if numel(firstLickPos) > 0
        firstLicks(t) = firstLickPos;
    else
        firstLicks(t) = NaN;
    end
end
clear('firstLickPos')
% [firstLicksSorted, indexes] = sort(firstLicks);
% shifts = abs(firstLicksSorted - lastTrackBin);
shifts = abs(firstLicks - lastTrackBin);

%licksInBinSorted = licksInBin(indexes,:);
%lickFreqInBinSorted = lickFreqInBin(indexes,:);
%nextHBStart = binNumber - homeBoxLength/binSize +1;
%binnedRoisExtDff = [binnedRoisDff nan(size(binnedRoisDff))];

binnedRoisLickAlignedDff = nan(size(binnedRoisDff).*[1 2 1]);
binnedRoisLickAlignedDeconv = nan(size(binnedRoisDeconv).*[1 2 1]);
binnedRoisLickAlignedDeconvRate = nan(size(binnedRoisDeconvRate).*[1 2 1]);


for roi = 1:1:nROIs
    
    for t = 1:1:allTrials
        if ~isnan(shifts(t))
           binnedRoisLickAlignedDff(t,1+shifts(t):binNumber+shifts(t),roi) = binnedRoisDff(t,1:binNumber,roi);
           binnedRoisLickAlignedDeconv(t,1+shifts(t):binNumber+shifts(t),roi) = binnedRoisDeconv(t,1:binNumber,roi);
           binnedRoisLickAlignedDeconvRate(t,1+shifts(t):binNumber+shifts(t),roi) = binnedRoisDeconvRate(t,1:binNumber,roi);
        end
    end
end

% DO THE SAME ALIGNMENT FOR BEHAVIOR TRIAL MATRICES

binVel = sData.behavior.trialMatrices.binVel;
licksInBin = sData.behavior.trialMatrices.licksInBin;
lickFreqInBin = sData.behavior.trialMatrices.lickFreqInBin;

binVelLickAligned = nan(size(binVel).*[1 2]);
licksInBinLickAligned = nan(size(licksInBin).*[1 2]);
lickFreqInBinLickAligned = nan(size(lickFreqInBin).*[1 2]);

for t = 1:1:allTrials
    if ~isnan(shifts(t))
        binVelLickAligned(t,1+shifts(t):binNumber+shifts(t)) = binVel(t,1:binNumber);
        licksInBinLickAligned(t,1+shifts(t):binNumber+shifts(t)) = licksInBin(t,1:binNumber);
        lickFreqInBinLickAligned(t,1+shifts(t):binNumber+shifts(t)) = lickFreqInBin(t,1:binNumber);
    end
end

%{
figure
hold on
c = zeros(size(licksInBinLickAligned));
c(licksInBinLickAligned > 0) = 1;
mymap = [1 1 1; 0 0 0];
imagesc(c)
colormap(gca,mymap)

figure
imagesc(binVelLickAligned)

figure
plot(binVelLickAligned')

figure
plot(nanmean(binVelLickAligned))
%}

sData.imdata(i).binnedRoisDff = binnedRoisDff;
sData.imdata(i).binnedRoisDeconv = binnedRoisDeconv;
sData.imdata(i).binnedRoisDeconvRate = binnedRoisDeconvRate;

sData.imdata(i).binnedRoisLickAlignedDff = binnedRoisLickAlignedDff;
sData.imdata(i).binnedRoisLickAlignedDeconv = binnedRoisLickAlignedDeconv;
sData.imdata(i).binnedRoisLickAlignedDeconvRate = binnedRoisLickAlignedDeconvRate;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create position tuning curves

binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
% binNumber = binNumber + homeBoxLength/binSize;

%%%%%% --- DFF ---
%%%% Trial Type 1
nTrialTypes = size(sData.trials.trialTypesMeta,2);

for t = 1:1:nTrialTypes
%avBinnedRoisTrialType1Dff(nROIs,binNumber) = NaN;
for roi = 1:1:nROIs 
    avBinnedRoisDff{t}(roi,:) = nanmean(binnedRoisDff(sData.trials.trialTypesMeta(t).trials,:,roi));
    avBinnedRoisDeconv{t}(roi,:) = nanmean(binnedRoisDeconv(sData.trials.trialTypesMeta(t).trials,:,roi));
    avBinnedRoisDeconvRate{t}(roi,:) = nanmean(binnedRoisDeconvRate(sData.trials.trialTypesMeta(t).trials,:,roi));
end

end

% Lick aligned
for t = 1:1:nTrialTypes
for roi = 1:1:nROIs 
    avBinnedRoisLickAlignedDff{t}(roi,:) = nanmean(binnedRoisLickAlignedDff(sData.trials.trialTypesMeta(t).trials,:,roi));
    avBinnedRoisLickAlignedDeconv{t}(roi,:) = nanmean(binnedRoisLickAlignedDeconv(sData.trials.trialTypesMeta(t).trials,:,roi));
    avBinnedRoisLickAlignedDeconvRate{t}(roi,:) = nanmean(binnedRoisLickAlignedDeconvRate(sData.trials.trialTypesMeta(t).trials,:,roi));
end

end

sData.imdata(i).avBinnedRois.avBinnedRoisDff = avBinnedRoisDff;
sData.imdata(i).avBinnedRois.avBinnedRoisDeconv = avBinnedRoisDeconv;
sData.imdata(i).avBinnedRois.avBinnedRoisDeconvRate = avBinnedRoisDeconvRate;
sData.imdata(i).avBinnedRois.avBinnedRoisLickAlignedDff = avBinnedRoisLickAlignedDff;
sData.imdata(i).avBinnedRois.avBinnedRoisLickAlignedDeconv = avBinnedRoisLickAlignedDeconv;
sData.imdata(i).avBinnedRois.avBinnedRoisLickAlignedDeconvRate = avBinnedRoisLickAlignedDeconvRate;

end

end