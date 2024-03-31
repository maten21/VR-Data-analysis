% Calculate spatial modulation Z-score Do shuffling on the entire 3D binned
% ROI array
% 
% Input: trial matrix: M(trials,positionBins,ROI)
%           signThreshold default = 0.995;
% Outputs: spatial modulation index (z-score), location of peak 
%
% M = sData.imdata.binnedRoisDeconvRate(sData.trials.contextsMeta(1).trials,:,48);
% M = sData.imdata.binnedRoisDff(sData.trials.contextsMeta(1).trials,:,38);

function [corrCoef, isSignificant] = getSMICorrAllRois(M1,M2,signThreshold)
%???? What should I do with the NaNs?
if nargin < 3
    signThreshold = 0.995;
end

smoothSpan = 5;

M1 = fillmissing(M1,'linear',2);
MRand1 = M1;
M2 = fillmissing(M2,'linear',2);
MRand2 = M2;


corrCoefShuf = nan(1000,1);

for j = 1:1:1000

% randomize matrix
for r = 1:1:size(M1,1)
    temp = M1(r,:,:);
    i = randi(size(M1,2));
    
    if i ~= 1
        MRand1(r,:,:) = [temp(1,i:size(M1,2),:) temp(1,1:i-1,:)];
    else
        MRand1(r,:,:) = temp;
    end
end
% randomize matrix
for r = 1:1:size(M2,1)
    temp = M2(r,:,:);
    i = randi(size(M2,2));
    
    if i ~= 1
        MRand2(r,:,:) = [temp(1,i:size(M2,2),:) temp(1,1:i-1,:)];
    else
        MRand2(r,:,:) = temp;
    end
end

% even = 3:2:size(M,1); 
% odd = 2:2:size(M,1); 

% Analyze tuning curve

tunCurves1 = smoothdata(permute(nanmean(MRand1),[2 3 1]),'gaussian',smoothSpan);
tunCurves2 = smoothdata(permute(nanmean(MRand2),[2 3 1]),'gaussian',smoothSpan);

for i = 1:1:size(tunCurves1,2)
    corrCoefShuf(j,i) = corr(tunCurves1(:,i),tunCurves2(:,i));
end

end

tunCurve1 = smoothdata(permute(nanmean(M1),[2 3 1]),'gaussian',smoothSpan);
tunCurve2 = smoothdata(permute(nanmean(M2),[2 3 1]),'gaussian',smoothSpan);

for i = 1:1:size(tunCurve1,2)
    
    corrCoef(i) = corr(tunCurve1(:,i),tunCurve2(:,i));
    isSignificant(i) = corrCoef(i) > quantile(corrCoefShuf(:,i),signThreshold);

end





end

%{

% even = 3:2:size(M,1); 
% odd = 2:2:size(M,1); 

M = smoothdata(sData.imdata.binnedRoisDff(sData.trials.contextsMeta(4).trials,:,77),2,'gaussian',5);
M= M(2:end,:);


figure
subplot(1,3,1)
imagesc(M(even,:))

subplot(1,3,2)
imagesc(M(odd,:))


subplot(1,3,3)
hold on
plot(nanmean(M(even,:)))
plot(nanmean(M(odd,:)))

corr(nanmean(M(even,:))',nanmean(M(odd,:))')



Mav = smoothdata(sData.imdata.binnedRoisDff(sData.trials.contextsMeta(1).trials,:,:),2,'gaussian',5);

Ma = nanmean(Mav,3);

figure
imagesc(Ma)

figure
plot(nanmean(Ma))

figure
histogram(corrCoefShuf)

quantile(corrCoefShuf,0.99)
quantile(corrCoefShuf,0.995)
%}

%% Plot figure

%{
figure; 

subplot(1,2,1); 
hold on
imagesc(M);
caxis([0 15])
plot(smoothdata(mean(M),2,'gaussian',smoothSpan)*15,'w-','LineWidth',3)
axis([0, size(M,2), 0, size(M,1)])

subplot(1,2,2);
hold on
imagesc(MRand); 
caxis([0 15])
plot(smoothdata(mean(MRand),2,'gaussian',smoothSpan)*15,'w-','LineWidth',3)
axis([0, size(M,2), 0, size(M,1)])

%}













%{
figure; histogram(peakAmplRand)
peakAmplRand

figure; imagesc(tunCurves)
figure; imagesc(corr(tunCurves))
colormap('jet')

figure; hold on; plot(tunCurve); plot(tunCurves(1:10,:)')
%}

%{
function [SMI, peakAmplRand, peakAmplExp] = getSMP(M,fraction)

M = binnedRoisDeconvRate(trials,:,129); 

M = fillmissing(M,'linear',2);
MRand = M;


smoothSpan = 9;
peakAmplRand = NaN(1000,1);
peakAmplExp = NaN(1000,1);
%tunCurves = NaN(1000,size(M,2));

for j = 1:1:100

% randomize matrix
for r = 1:1:size(M,1)
    temp = M(r,:);
    i = randi(size(M,2));
    
    if i ~= 1
        MRand(r,:) = [temp(i:size(M,2)) temp(1:i-1)];
    else
        MRand(r,:) = temp;
    end
end

% Analyze tuning curve

for k = 1:1:10
tr = randperm(size(M,1),floor(size(M,1)*fraction));

tunCurveR = smoothdata(nanmean(MRand(tr,:)),2,'gaussian',smoothSpan);
peakAmplRand((j-1)*10+k) = (max(tunCurveR) - min(tunCurveR)) / mean(tunCurveR);

tunCurveE = smoothdata(nanmean(M(tr,:)),2,'gaussian',smoothSpan);
peakAmplExp((j-1)*10+k) = (max(tunCurveE) - min(tunCurveE)) / mean(tunCurveE);

end

end


figure; hold on; histogram(peakAmplRand); histogram(peakAmplExp);

%peakBin = find(tunCurveE == max(tunCurveE));

SMI = %% abs(nanmean(peakAmplExp)-nanmean(peakAmplRand))/std(peakAmplRand);

end
%}


