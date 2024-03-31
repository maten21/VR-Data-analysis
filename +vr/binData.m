

function [M, Xax] = binData(sData,dataSignal,binSize,binSignal)

if nargin < 4 % position is default, could be time, etc.
    if length(dataSignal) == length(sData.behavior.signals.corridorPositionDs)
        binSignal = sData.behavior.signals.corridorPositionDs;        
    else length(dataSignal) == length(sData.behavior.signals.corridorPosition)
        binSignal = sData.behavior.signals.corridorPosition; 
    end
end

if nargin < 3
    binSize = 2;
end

if length(dataSignal) == length(sData.behavior.signals(1).trialIndexesDs)
    trialIndexes = sData.behavior.signals(1).trialIndexesDs;
else length(dataSignal) == length(sData.behavior.signals(1).trialIndexes)
    trialIndexes = sData.behavior.signals(1).trialIndexes;
end
    



minBin = floor(min(binSignal));
maxBin = ceil(max(binSignal));
rangeBin = maxBin-minBin;


trialStartIndexes = find(diff(binSignal) < -rangeBin/10) +1; 
nAllTrials = numel(trialStartIndexes)-1;

Xax = (minBin+binSize) : binSize : maxBin;

pos = discretize(binSignal,Xax);


M = nan(nAllTrials,rangeBin/binSize);
for i = 1:1:nAllTrials

s = dataSignal(trialIndexes==i);
p = pos(trialIndexes==i);

for j = 1:1:rangeBin/binSize
M(i,j) = nanmean(s(p==j));
end

end




end