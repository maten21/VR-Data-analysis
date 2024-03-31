function [] = calcSessionMedianSfifts(fileName,filePath)
% This function calculates the shift of the velocity and lick curves
% of different trial tyles compared to control (trial type 1)
%
%
%

load(fullfile(filePath,fileName));

binSize = sData.behavior.trialMatrices.meta.binSize;
firstBin = find(sData.stats.sessionAvs(1).plotXAxis == binSize);    %60/2 + 1; % first bin of the corridor
lastBin = firstBin + sData.behavior.rewardZone/binSize - 1; % last bin before RZ
centerBin = ceil(mean([firstBin lastBin])); % central bin

for t = 1:1:size(sData.trials.trialTypesMeta,2)

    velocitys(t,:) = sData.stats.sessionMedians(t).medBinVel;
    lickFreqs(t,:) = sData.stats.sessionMedians(t).medLickFreqInBin;
    
end
    
for i = 1:1:21 % generate matrix from shifted control curves
    
    shiftVels(i,i:numel(velocitys(1,:))+i-1) = velocitys(1,:);
    shiftLickFreqs(i,i:numel(lickFreqs(1,:))+i-1) = lickFreqs(1,:);
    %shiftVels(i,i:numel(velocitys(1,centerBin:lastBin))+i-1) = velocitys(1,centerBin:lastBin);
    %shiftLickFreqs(i,i:numel(lickFreqs(1,centerBin:lastBin))+i-1) = lickFreqs(1,centerBin:lastBin);
    
end

shift = -20:2:20;


for i = 1:1:51 % generate matrix from shifted control curves
    
    shiftVels(i,i:numel(velocitys(1,:))+i-1) = velocitys(1,:);
    shiftLickFreqs(i,i:numel(lickFreqs(1,:))+i-1) = lickFreqs(1,:);
    %shiftVels(i,i:numel(velocitys(1,centerBin:lastBin))+i-1) = velocitys(1,centerBin:lastBin);
    %shiftLickFreqs(i,i:numel(lickFreqs(1,centerBin:lastBin))+i-1) = lickFreqs(1,centerBin:lastBin);
    
end

shift = -50:2:50;



for t = 1:1:size(sData.trials.trialTypesMeta,2)
    
    corrVel = corr(shiftVels(:,centerBin+10:lastBin+10)',velocitys(t,centerBin:lastBin)');
    corrFreq = corr(shiftLickFreqs(:,centerBin+10:lastBin+10)',lickFreqs(t,centerBin:lastBin)');
    
    sData.trials.trialTypesMeta(t).velShift = shift(find(corrVel == max(corrVel)));
    sData.trials.trialTypesMeta(t).lickFreqShift = shift(find(corrFreq == max(corrFreq)));    
    
end

save(fullfile(filePath,fileName),'sData');
clear('sData');

end


Z = [0 0 0 0 0 0 0 0 0 0];
longFreq = [Z lickFreqs(t,:) Z];
longVel = [Z velocitys(t,:) Z];



    corrVel = corr(shiftVels(:,11:size(shiftVels,2)-10)',velocitys(t,centerBin:lastBin)');
    corrFreq = corr(shiftLickFreqs(:,centerBin+10:lastBin+10)',lickFreqs(t,centerBin:lastBin)');

corrVel = corr(shiftVels(:,centerBin+25:lastBin+25)',velocitys(t,centerBin:lastBin)','Type','Spearman');
corrFreq = corr(shiftLickFreqs(:,centerBin+25:lastBin+25)',lickFreqs(t,centerBin:lastBin)','Type','Spearman');

corr([1 2 3]',[NaN 4 5]')

figure
hold on
plot(shiftVels(:,centerBin+25:lastBin+25)')
plot(velocitys(t,centerBin:lastBin)','k')
plot(shiftVels(23,centerBin+25:lastBin+25)','r')

figure
plot(velocitys(t,centerBin:lastBin)',shiftVels(23,centerBin+25:lastBin+25)','o')

figure
hold on
plot(velocitys(t,centerBin:lastBin)'-shiftVels(23,centerBin+25:lastBin+25)')
plot(velocitys(t,centerBin:lastBin)')
plot(shiftVels(23,centerBin+25:lastBin+25)')







y = resample(x,ups,dns);

V1 = velocitys(1,1:lastBin);
V2 = velocitys(2,1:lastBin);

L1 = lickFreqs(1,1:lastBin);
L2 = lickFreqs(2,1:lastBin);



j = 1;
for i = 70:1:110
nBins(j) = numel(resample(V1,i,1));
    V1Stretch(j,1:nBins(j)) = resample(V1,i,1);
    L1Stretch(j,1:nBins(j)) = resample(L1,i,1);
    j=j+1;    
end
stretchFactor = 70:1:110;
V2Stretch = resample(V2,100,1);
L2Stretch = resample(L2,100,1);

corrAllBins = corr(L1Stretch(:,1:numel(L2Stretch))',L2Stretch');
corrAllBinsS = corr(L1Stretch(:,1:numel(L2Stretch))',L2Stretch','Type','Spearman');


figure
hold on
plot(stretchFactor,corrAllBins)
plot(stretchFactor,corrAllBinsS)


% THIS IS THE ACCURATE CALCULATION
for i = 1:1:numel(stretchFactor)
    shorter = min([nBins(i) numel(V2Stretch)]);
    corrVel(i) = corr(V1Stretch(i,1:shorter)',V2Stretch(1:shorter)');    
    corrVelS(i) = corr(V1Stretch(i,1:shorter)',V2Stretch(1:shorter)','Type','Spearman');  
    corrLick(i) = corr(L1Stretch(i,1:shorter)',L2Stretch(1:shorter)');    
    corrLickS(i) = corr(L1Stretch(i,1:shorter)',L2Stretch(1:shorter)','Type','Spearman');  
end

stretchFactor(find(corrVel == max(corrVel)));
stretchFactor(find(corrVelS == max(corrVelS)));
stretchFactor(find(corrLick == max(corrLick)));
stretchFactor(find(corrLickS == max(corrLickS)));


figure % corr values
hold on
plot(stretchFactor,corrLick)
plot(stretchFactor,corrLickS)


figure
hold on
plot(1:numel(L1Stretch(31,:)),L1Stretch(31,:))
plot(1:numel(L2Stretch),L2Stretch)
plot(1:numel(L1Stretch(find(corrShorter == max(corrShorter)),:)),L1Stretch(find(corrShorter == max(corrShorter)),:))
%plot(1:numel(V1Stretch(find(corrShorterS == max(corrShorterS)),:)),V1Stretch(find(corrShorterS == max(corrShorterS)),:))



figure
plot(L1Stretch(i,1:nBins(17))',L2Stretch(1:nBins(17))','.')
plot(L1Stretch(find(corrShorter == max(corrShorter)),:)')
