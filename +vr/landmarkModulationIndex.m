

%% Calculate landmark modulation index (LMI)

function sData = landmarkModulationIndex(sData)

for t = 1:1:max(size(sData.trials.trialTypesMeta))
       
    trials = sData.trials.trialTypesMeta(t).trials;
    binSize = sData.behavior.trialMatrices.meta.binSize;
    firstLickPositionCm = sData.stats.lickDistribution.firstLickPositionCm';
    
    M = sData.behavior.trialMatrices.lickFreqInBin; % for calculating start positions
    firstBins = nan(size(M(:,1))); % for calculating start positions
        
    for r = 1:1:size(M,1) % calculating start positions        
        firstBins(r) = min(find(~isnan(M(r,:))));       
    end
    
    discr = discretize([firstBins', 1, 31],4);
    category = discr(1:numel(discr)-2)';
    

    randStartPos(4) = struct();
    
        firstBinsCm = (firstBins-1)*binSize -140;
        firstLickFromStartCm = firstLickPositionCm - firstBinsCm;
        % firstLickFromStartCm = firstLickPositionCm + 140 - (firstBins-1)*binSize;
        
    for i = 1:1:4
        
        trialSubType = intersect(trials,find(category==i));
        
        % Save these data to sData
        randStartPos(i).trials = trialSubType;
        randStartPos(i).firstLickPositionCm = nanmedian(firstLickPositionCm(trialSubType));
        randStartPos(i).firstLickFromStartCm = nanmedian(firstLickFromStartCm(trialSubType));
        randStartPos(i).startPosCm = nanmean(firstBinsCm(trialSubType));
        
        lickPositions(i) = randStartPos(i).firstLickPositionCm;
        licksFromStart(i) = randStartPos(i).firstLickFromStartCm;
        startPosCm(i) = randStartPos(i).startPosCm;
        
    end
    
    
    % LMI = std(lickPositions)/std(licksFromStart)
    %LMI = std(licksFromStart)/std(startPosCm);
    
    sData.trials.trialTypesMeta(t).randStartPos = randStartPos;
    
    
    clear('randStartPos')
        %sData.trials.trialTypesMeta(t).firstLicksPerc = quantile(sData.stats.lickDistribution.firstLickPositionPerc(trials),[0.25 0.50 0.75]);
    
        
        for j = 1:1:1000
            
            randCategory = randi(4,numel(category)); % randomize
            
                
            for i = 1:1:4
                
                trialSubType = intersect(trials,find(randCategory==i));
                % randStartPos(i).trials = trialSubType;
                randFirstLickPositionCm(j,i) = nanmedian(firstLickPositionCm(trialSubType));
                randFirstLickFromStartCm(j,i) = nanmedian(firstLickFromStartCm(trialSubType));
                %randStartPos(i).startPosCm = nanmean(firstBinsCm(trialSubType));
                
                %lickPositions(i) = randStartPos(i).firstLickPositionCm;
                %licksFromStart(i) = randStartPos(i).firstLickFromStartCm;
                %startPosCm(i) = randStartPos(i).startPosCm;
                
            end
   
        end


        landmarkModZscore(t) = abs(nanstd(licksFromStart) - mean(nanstd(randFirstLickFromStartCm,0,2))) / std(nanstd(randFirstLickFromStartCm,0,2));
        startModZscore(t) = abs(nanstd(lickPositions) - mean(nanstd(randFirstLickPositionCm,0,2))) / std(nanstd(randFirstLickPositionCm,0,2));
        
        LMI = landmarkModZscore(t)./(landmarkModZscore(t) + startModZscore(t));
        
        sData.trials.trialTypesMeta(t).landModInd = LMI;
end

end

%{

figure
subplot(1,2,1)
hold on
histogram(nanstd(randFirstLickFromStartCm,0,2))
histogram(nanstd(licksFromStart))

subplot(1,2,2)
hold on
histogram(nanstd(randFirstLickPositionCm,0,2))
histogram(nanstd(lickPositions))


landmarkModZscore = abs(nanstd(licksFromStart) - mean(nanstd(randFirstLickFromStartCm,0,2))) / std(nanstd(randFirstLickFromStartCm,0,2))
startModZscore = abs(nanstd(lickPositions) - mean(nanstd(randFirstLickPositionCm,0,2))) / std(nanstd(randFirstLickPositionCm,0,2))


LMI
LMIZ = landmarkModZscore./(landmarkModZscore + startModZscore)


CovStart = std(nanstd(randFirstLickFromStartCm,0,2)) / mean(nanstd(randFirstLickFromStartCm,0,2))
CovLandm = std(nanstd(randFirstLickPositionCm,0,2)) / mean(nanstd(randFirstLickPositionCm,0,2))



figure; 


errorbar(nanmean(firstLickFromStartCm),nanstd(firstLickFromStartCm),'o')


%}

