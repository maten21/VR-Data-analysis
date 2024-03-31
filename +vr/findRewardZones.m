
function rewardZones = findRewardZones(sData)

binSize = sData.behavior.trialMatrices.meta.binSize;
rewardZones = [];


for context = 1:1:length(sData.trials.contextsMeta)
    rew = mean(sData.behavior.trialMatrices.rewardInBin(sData.trials.contextsMeta(context).trials,:));
    rew(isnan(rew)) = 0;
    rewardZ = sData.stats.sessionAvs(1).plotXAxis(find(rew));
    
    rewardZones(context,1) = rewardZ(1);
    rewardZones(context,2) = rewardZ(find(diff(discretize(rewardZ,2)))+1);
    
    
end
rewardZones = rewardZones - binSize;


end

