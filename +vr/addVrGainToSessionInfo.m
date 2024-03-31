% Check the gain of the recording

function addVrGainToSessionInfo(fileName,filePath)

load(fullfile(filePath,fileName));


trialStartIndexes = find(diff(sData.behavior.signals.corridorPosition) < -100) +1;

gains(1:1:numel(trialStartIndexes)-1) = NaN;

for i = 1:numel(trialStartIndexes)-1
    j = trialStartIndexes(i);
    k = trialStartIndexes(i+1)-1;
    
    vrDistance = sData.behavior.signals.corridorPosition(k) - sData.behavior.signals.corridorPosition(j);
    wheelDistance = sData.daqdata.wheelDistance(k) - sData.daqdata.wheelDistance(j);
    gains(i) = vrDistance/wheelDistance;
    
end

sData.sessionInfo.vrGain = nanmean(gains);

save(fullfile(filePath,fileName),'sData');
clear('sData');

end





