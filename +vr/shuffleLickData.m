%% bootstrapping



% clear
function [hitRateIndexZ, hitRateRand, hitRateExp, spatialModIndZ, SMIrand , SMIexp ] = shuffleLickData(sData,trialSubset,nTrialsIncl)

binSize = 2;
plotXAxis = sData.stats.sessionAvs.plotXAxis;

licksInBin = sData.behavior.trialMatrices.licksInBin(trialSubset,:);
rewardInBin = sData.behavior.trialMatrices.rewardInBin(trialSubset,:);
nRows = size(licksInBin,1);

minPos = find(plotXAxis == binSize);
maxPos = find(plotXAxis == max(plotXAxis)); % CHECK THIS LATER!!! Exclude drinking!

for r = 1:1:nRows  % Exclude drinking!
    rewGiven = false;
    for c = maxPos-4:1:maxPos
        if rewGiven
            licksInBin(r,c) = 0;
        %elseif ~rewGiven && rewardInBin(r,c) > 0 && licksInBin(r,c) > 1
        %    licksInBin(r,c) = 1;
        %    rewGiven = true;
        elseif ~rewGiven && rewardInBin(r,c) > 0 
            licksInBin(r,c) = 1;
            rewGiven = true;
        end
    end
end

%{
for i = 1:1:nRows
    posMatrixFull(i,:) = plotXAxis;

end
%}


%lickMatrix = licksInBin(:,j:k);
%posMatrix = posMatrixFull(:,j:k);

maxPosIndex = max(plotXAxis)/binSize;


hitRateRand(1:1000,1) = NaN;
hitRateExp(1:1000,1) = NaN;
SMIrand(1:1000,1) = NaN;
SMIexp(1:1000,1) = NaN;
% original hit rate 


%tic;
% Randomize
for r = 1:1:1000

% this is to tormalize for the trial number
randSubset = randperm(numel(trialSubset),nTrialsIncl);    
lickMatrix = licksInBin(randSubset,minPos:maxPos);
   
nRows = nTrialsIncl;
    
R = randi(maxPosIndex-1,[nRows,1]);
%R = reshape(R,[10 10]);

licksLin = lickMatrix';
licksLin = licksLin(:);

indexes(1:nRows*2,1:2) = NaN;
for i = 2:2:nRows*2
    j = i/2;
    indexes(i-1,1) = 1 + j*maxPosIndex-maxPosIndex;
    indexes(i-1,2) = R(j) + j*maxPosIndex-maxPosIndex;
    indexes(i,1) = R(j)+1 + j*maxPosIndex-maxPosIndex;
    indexes(i,2) = maxPosIndex + j*maxPosIndex-maxPosIndex;
end

randRows = randperm(nRows*2);
licksLinRand = licksLin;
%licksLinRand = NaN;

m = 1;
for i = 1:1:(nRows*2)
    j = indexes(randRows(i),1);
    k = indexes(randRows(i),2);
    temp = licksLin(j:k);
    n = m + numel(temp) - 1;
    licksLinRand(m:n) = temp;
    m = n + 1;
    
end

licksLinRand = licksLinRand';

lickMatrixRand = reshape(licksLinRand,[numel(licksLinRand)/nRows nRows]);
lickMatrixRand = lickMatrixRand';
    
hitRateRand(r) = sum(sum(lickMatrixRand(:,maxPosIndex-4:maxPosIndex),2)>0)/nRows;

tuningCurveRand = smoothdata(sum(lickMatrixRand),'gaussian',9);

SMIrand(r) = (max(tuningCurveRand) - min(tuningCurveRand))/mean(tuningCurveRand);


hitRateExp(r) = sum(sum(lickMatrix(:,maxPosIndex-4:maxPosIndex),2)>0)/nRows;
tuningCurveExp = smoothdata(sum(lickMatrix),'gaussian',9);
SMIexp(r) = (max(tuningCurveExp) - min(tuningCurveExp))/mean(tuningCurveExp);    
    

end
%toc;



hitRateIndexZ = (mean(hitRateExp)-mean(hitRateRand))/std(hitRateRand); % This is the Z score!
spatialModIndZ = (mean(SMIexp)-mean(SMIrand))/std(SMIrand); % This is the Z score!
%SMIs = hitRate./hitRateRand;

%SMI = hitRate/mean(hitRateRand);
%SMIs = hitRate./hitRateRand;

end

%{

%figure; histogram(hitRateExp)
figure; 
hold on
h = histogram(SMIrand);
ax = gca;
rectangle('Position',[SMIexp-h.BinWidth/4, ax.YLim(1), h.BinWidth/2, ax.YLim(2)],'FaceColor',[1 0.1 0.1],'EdgeColor','none'); % home-box



figure
subplot(1,2,1)

imagesc(lickMatrix*-1+1)
colormap('gray')
caxis([0 1])


subplot(1,2,2)

imagesc(lickMatrixRand*-1+1)
colormap('gray')
caxis([0 1])





figure
subplot(1,2,1)

plot(smoothdata(sum(lickMatrix),'gaussian',9))
%colormap('gray')
%caxis([0 1])


subplot(1,2,2)

plot(smoothdata(sum(lickMatrixRand),'gaussian',9))
%colormap('gray')
%caxis([0 1])
%}