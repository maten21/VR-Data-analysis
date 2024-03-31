

function bootStrapData = isSignificantBootstrap(trialMatrix,optoTrials)

shufflingFactor = 10000;


nOptoTrials = length(optoTrials);
nAllTrials = size(trialMatrix,1);
ctrlTrials = setdiff(1:1:nAllTrials,optoTrials); 

ampl1 = nan(shufflingFactor,1);
ampl2 = nan(shufflingFactor,1);

for i = 1:1:shufflingFactor

randTrials2 = randperm(nAllTrials,nOptoTrials); 
randTrials1 = setdiff(1:1:nAllTrials,randTrials2); 

tuneCurve1 = nanmedian(trialMatrix(randTrials1,:));
tuneCurve2 = nanmedian(trialMatrix(randTrials2,:));

% Could be refined with better baseline and amplitude method
ampl1(i) = max(tuneCurve1) - min(tuneCurve1);
ampl2(i) = max(tuneCurve2) - min(tuneCurve2);



end


tuneCurveCtrl = nanmedian(trialMatrix(ctrlTrials,:));
tuneCurveOpto = nanmedian(trialMatrix(optoTrials,:));

amplCtrl = max(tuneCurveCtrl) - min(tuneCurveCtrl);
amplOpto = max(tuneCurveOpto) - min(tuneCurveOpto);


effectRand = ampl2 - ampl1;
effect = amplOpto - amplCtrl;

p1 = sum(effect < effectRand)/length(effectRand);
p2 = sum(effect > effectRand)/length(effectRand);

bootStrapData.amplCtrl = amplCtrl;
bootStrapData.amplOpto = amplOpto;
bootStrapData.ampl1 = ampl1;
bootStrapData.ampl2 = ampl2;
bootStrapData.p = min(p1,p2);


end

