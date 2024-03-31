function [] = calcSessionMedianStretch(fileName,filePath)
% This function calculates the shift of the velocity and lick curves
% of different trial tyles compared to control (trial type 1)
%
%
%

load(fullfile(filePath,fileName));

binSize = sData.behavior.trialMatrices.meta.binSize;
firstBin = find(sData.stats.sessionAvs(1).plotXAxis == binSize);    %60/2 + 1; % first bin of the corridor
rewardZone = sData.behavior.rewardZone;
lastBin = firstBin + rewardZone/binSize - 1; % last bin before RZ
%centerBin = ceil(mean([firstBin lastBin])); % central bin


for t = 1:1:size(sData.trials.trialTypesMeta,2)

    ctrlVel = sData.stats.sessionMedians(1).medBinVel(1:lastBin);
    ctrlLickFreq = sData.stats.sessionMedians(1).medLickFreqInBin(1:lastBin);
    vel = sData.stats.sessionMedians(t).medBinVel(1:lastBin)';
    lickFreq = sData.stats.sessionMedians(t).medLickFreqInBin(1:lastBin)';
    
% calculate stretched curves
j = 1;
for i = 70:1:110
lastBins(j) = numel(resample(ctrlVel,i,1));
    stretchedCtrlVel(j,1:lastBins(j)) = resample(ctrlVel,i,1);
    stretchedCtrlLickFreq(j,1:lastBins(j)) = resample(ctrlLickFreq,i,1);
    j = j + 1;    
end

vel = resample(vel,100,1);
lickFreq = resample(lickFreq,100,1);

stretchFactor = 70:1:110;


% THIS IS THE ACCURATE CALCULATION
for i = 1:1:numel(stretchFactor)
    shorter = min([lastBins(i) numel(vel)]);
    corrVel(i) = corr(stretchedCtrlVel(i,1:shorter)',vel(1:shorter)');    
    %corrVelS(i) = corr(stretchedCtrlVel(i,1:shorter)',vel(1:shorter)','Type','Spearman');  
    corrLick(i) = corr(stretchedCtrlLickFreq(i,1:shorter)',lickFreq(1:shorter)');    
    %corrLickS(i) = corr(stretchedCtrlLickFreq(i,1:shorter)',lickFreq(1:shorter)','Type','Spearman');  
end

stretchFactor(find(corrVel == max(corrVel)));
%stretchFactor(find(corrVelS == max(corrVelS)));
stretchFactor(find(corrVel == max(corrVel)));
%stretchFactor(find(corrLickS == max(corrLickS)));


sData.trials.trialTypesMeta(t).velStretchFactor = stretchFactor(find(corrVel == max(corrVel)));
sData.trials.trialTypesMeta(t).lickStretchFactor = stretchFactor(find(corrLick == max(corrLick)));  

sData.trials.trialTypesMeta(t).velShift = (stretchFactor(find(corrVel == max(corrVel)))-100) * rewardZone/100;
sData.trials.trialTypesMeta(t).lickFreqShift = (stretchFactor(find(corrLick == max(corrLick)))-100) * rewardZone/100;   
    





end



save(fullfile(filePath,fileName),'sData');
clear('sData');


end