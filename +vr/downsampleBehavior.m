 function sData = downsampleBehavior(sData)

try
    
corridorPosition = sData.behavior.signals(1).corridorPosition; 
trialStartIndexes = find(diff(corridorPosition) < -10)+1;
nAllTrials = numel(trialStartIndexes)-1;

lickPositions = sData.behavior.signals(1).lickPositions;

fs = sData.daqdata.meta.fs;

positionRelativeToStart = nan(size(lickPositions));
positionRelativeToFirstLick = nan(size(lickPositions));
timeRelativeToStart = nan(size(lickPositions));
timeRelativeToLandmark = nan(size(lickPositions));
timeRelativeToFirstLick = nan(size(lickPositions));


for i = 1:1:nAllTrials
    
    j = trialStartIndexes(i);
    k = trialStartIndexes(i+1)-1;
    
    tempPos = corridorPosition(j:k);
    tempLickPos = lickPositions(j:k);
    temp = tempLickPos(tempLickPos>0);
        
    tempTime = linspace(0,numel(tempPos)/fs,numel(tempPos));
    landmarkIndex = find(abs(tempPos) == min(abs(tempPos)));
    landmarkIndex = landmarkIndex(1);
    timeRelativeToStart(j:k) = tempTime;
    timeRelativeToLandmark(j:k) = tempTime - landmarkIndex/fs;
    positionRelativeToStart(j:k) = tempPos - tempPos(1);
    
    if numel(temp) > 0
        positionRelativeToFirstLick(j:k) = tempPos - temp(1);
        tempFirstLickIndex = find(tempLickPos == temp(1));
        timeRelativeToFirstLick(j:k) = tempTime - (tempFirstLickIndex(1)/fs);
    end
    
    clear('tempPos', 'tempLickPos', 'temp', 'tempTime', 'landmarkIndex');
end
    
catch
end


%% 

corridorPosition = sData.behavior.signals(1).corridorPosition; 

trialStartIndexes = find(diff(corridorPosition) < -10) +1; % changed from -100 to capture very first trials
trialIndexes = nan(size(sData.behavior.signals(1).corridorPosition));

fs = sData.daqdata.meta.fs;


for t = 1:1:numel(trialStartIndexes)-1
    trialIndexes(trialStartIndexes(t):trialStartIndexes(t+1)-1) = t;   
end

sData.behavior.signals(1).trialIndexes = trialIndexes;

%corridorPosition = sData.behavior.signals(1).corridorPosition; 



nFOVs =  length(sData.imdata);

P1FOVs = 0;
P2FOVs = 0;
for i = 1:1:nFOVs
    try
        if strcmp(sData.imdata(i).meta.imagingSystem(1:5),'Path1')
            P1FOVs = P1FOVs + 1;
        end
        if strcmp(sData.imdata(i).meta.imagingSystem(1:5),'Path2')
            P2FOVs = P2FOVs + 1;
        end
    catch
    end
end

P1Count = -1;
P2Count = -1;

for i = 1:1:nFOVs

try
    if isequal(sData.imdata(i).meta.imagingSystem(1:5),'Path1')
        scanFrameRate = sData.imdata(i).meta.scanFrameRate;
        samplePerFrame = fs/scanFrameRate;
        frameIndexes = vr.fixFrameIndexes(sData.daqdata.frameIndex, samplePerFrame);
    elseif isequal(sData.imdata(i).meta.imagingSystem(1:5),'Path2')
        scanFrameRate = sData.imdata(i).meta.scanFrameRate;
        samplePerFrame = fs/scanFrameRate;
        frameIndexes = vr.fixFrameIndexes(sData.daqdata.frameIndex2, samplePerFrame);
    end
catch
    scanFrameRate = sData.imdata.meta.fps;
    samplePerFrame = fs/scanFrameRate;
    frameIndexes = vr.fixFrameIndexes(sData.daqdata.frameIndex, samplePerFrame);
end



if numel(frameIndexes) > sData.imdata(i).nSamples
    frameIndexes = frameIndexes(1:sData.imdata(i).nSamples);
end

deltaT = mean(diff(frameIndexes))/fs;

try
    if strcmp(sData.imdata(i).meta.imagingSystem(1:5),'Path1')
        P1Count = P1Count + 1;
        timeDs = frameIndexes/fs + deltaT/P1FOVs*P1Count;
    end
    if strcmp(sData.imdata(i).meta.imagingSystem(1:5),'Path2')
        P2Count = P2Count + 1;
        timeDs = frameIndexes/fs + deltaT/P2FOVs*P2Count;
    end
catch
    timeDs = frameIndexes/fs;
end



%{
if isfield(sData.imdata(i).meta,'scanFrameRate')
    deltaTime = 1/sData.imdata(i).meta.scanFrameRate;
elseif isfield(sData.imdata(i).meta,'fps')
    deltaTime = 1/sData.imdata(i).meta.fps;
end
%}

sData.behavior.signals(i).trialIndexesDs = trialIndexes(frameIndexes);
sData.behavior.signals(i).corridorPositionDs = corridorPosition(frameIndexes);
%sData.behavior.signals.lickPositionsDs = sData.behavior.signals.lickPositions(frameIndexes);
%sData.behavior.signals.lickEventsDs = sData.behavior.signals.lickEvents(frameIndexes);
%sData.behavior.signals.rewardEventsDs = sData.behavior.signals.rewardEvents(frameIndexes);
sData.behavior.signals(i).velocityDs = sData.behavior.signals(1).velocity(frameIndexes);
sData.behavior.signals(i).lickIfreqDs = sData.behavior.signals(1).lickIfreq(frameIndexes);
try
sData.behavior.signals(i).stimProtsDs = sData.daqdata.stimProtIndicators(frameIndexes); 
catch
end
sData.behavior.signals(i).timeDs = timeDs; %0:deltaTime:deltaTime*(numel(frameIndexes)-1); %the most accutate way is to calculate from the fram signals


try
% sData.behavior.signals.timeRelativeToStart = timeRelativeToStart;
% sData.behavior.signals.timeRelativeToLandmark = timeRelativeToLandmark;
% sData.behavior.signals.timeRelativeToFirstLick = timeRelativeToFirstLick;
% sData.behavior.signals.positionRelativeToStart = positionRelativeToStart;
% sData.behavior.signals.positionRelativeToFirstLick = positionRelativeToFirstLick;

sData.behavior.signals(i).timeRelativeToStartDs = timeRelativeToStart(frameIndexes);
sData.behavior.signals(i).timeRelativeToLandmarkDs = timeRelativeToLandmark(frameIndexes);
sData.behavior.signals(i).timeRelativeToFirstLickDs = timeRelativeToFirstLick(frameIndexes);
sData.behavior.signals(i).positionRelativeToStartDs = positionRelativeToStart(frameIndexes);
sData.behavior.signals(i).positionRelativeToFirstLickDs = positionRelativeToFirstLick(frameIndexes);
sData.behavior.signals(i).rewardDs = sData.daqdata.waterValve(frameIndexes);
catch
end

end

% if isfield(sData.daqdata,'optoSignal')
%    sData.behavior.signals.optoSignalDs = sData.daqdata.optoSignal(frameIndexes);    
% end
    

% BACKUP before multi FOV upgrade
%{

try
trialStartIndexes = find(diff(sData.daqdata.unityPosition) < -10)+1;
nAllTrials = numel(trialStartIndexes)-1;

lickPositions = sData.behavior.signals.lickPositions;
corridorPosition = sData.behavior.signals.corridorPosition; 
fs = sData.daqdata.meta.fs;

positionRelativeToStart = nan(size(lickPositions));
positionRelativeToFirstLick = nan(size(lickPositions));
timeRelativeToStart = nan(size(lickPositions));
timeRelativeToLandmark = nan(size(lickPositions));
timeRelativeToFirstLick = nan(size(lickPositions));


for i = 1:1:nAllTrials
    
    j = trialStartIndexes(i);
    k = trialStartIndexes(i+1)-1;
    
    tempPos = corridorPosition(j:k);
    tempLickPos = lickPositions(j:k);
    temp = tempLickPos(tempLickPos>0);
        
    tempTime = linspace(0,numel(tempPos)/fs,numel(tempPos));
    landmarkIndex = find(abs(tempPos) == min(abs(tempPos)));
    landmarkIndex = landmarkIndex(1);
    timeRelativeToStart(j:k) = tempTime;
    timeRelativeToLandmark(j:k) = tempTime - landmarkIndex/fs;
    positionRelativeToStart(j:k) = tempPos - tempPos(1);
    
    if numel(temp) > 0
        positionRelativeToFirstLick(j:k) = tempPos - temp(1);
        tempFirstLickIndex = find(tempLickPos == temp(1));
        timeRelativeToFirstLick(j:k) = tempTime - (tempFirstLickIndex(1)/fs);
    end
    
    clear('tempPos', 'tempLickPos', 'temp', 'tempTime', 'landmarkIndex');
end
    
catch
end


frameIndexes = find(diff(sData.daqdata.frameIndex)==1);
if numel(frameIndexes) > sData.imdata.nSamples
    frameIndexes = frameIndexes(1:sData.imdata.nSamples);
end

deltaTime = 1/sData.imdata.meta.fps;
trialStartIndexes = find(diff(sData.daqdata.unityPosition) < -10) +1; % changed from -100 to capture very first trials

trialIndexes = nan(size(sData.behavior.signals.corridorPosition));


for t = 1:1:numel(trialStartIndexes)-1
    trialIndexes(trialStartIndexes(t):trialStartIndexes(t+1)-1) = t;   
end

sData.behavior.signals.trialIndexes = trialIndexes;
sData.behavior.signals.trialIndexesDs = trialIndexes(frameIndexes);
sData.behavior.signals.corridorPositionDs = sData.behavior.signals.corridorPosition(frameIndexes);
%sData.behavior.signals.lickPositionsDs = sData.behavior.signals.lickPositions(frameIndexes);
%sData.behavior.signals.lickEventsDs = sData.behavior.signals.lickEvents(frameIndexes);
%sData.behavior.signals.rewardEventsDs = sData.behavior.signals.rewardEvents(frameIndexes);
sData.behavior.signals.velocityDs = sData.behavior.signals.velocity(frameIndexes);
sData.behavior.signals.lickIfreqDs = sData.behavior.signals.lickIfreq(frameIndexes);
sData.behavior.signals.stimProtsDs = sData.daqdata.stimProtIndicators(frameIndexes); 
sData.behavior.signals.time = 0:deltaTime:deltaTime*(numel(frameIndexes)-1);


try
% sData.behavior.signals.timeRelativeToStart = timeRelativeToStart;
% sData.behavior.signals.timeRelativeToLandmark = timeRelativeToLandmark;
% sData.behavior.signals.timeRelativeToFirstLick = timeRelativeToFirstLick;
% sData.behavior.signals.positionRelativeToStart = positionRelativeToStart;
% sData.behavior.signals.positionRelativeToFirstLick = positionRelativeToFirstLick;

sData.behavior.signals.timeRelativeToStartDs = timeRelativeToStart(frameIndexes);
sData.behavior.signals.timeRelativeToLandmarkDs = timeRelativeToLandmark(frameIndexes);
sData.behavior.signals.timeRelativeToFirstLickDs = timeRelativeToFirstLick(frameIndexes);
sData.behavior.signals.positionRelativeToStartDs = positionRelativeToStart(frameIndexes);
sData.behavior.signals.positionRelativeToFirstLickDs = positionRelativeToFirstLick(frameIndexes);
sData.behavior.signals.rewardDs = sData.daqdata.waterValve(frameIndexes);
catch
end

% if isfield(sData.daqdata,'optoSignal')
%    sData.behavior.signals.optoSignalDs = sData.daqdata.optoSignal(frameIndexes);    
% end
    


%}



    
end