%% SPLIT SORTED ROIs TO SUBSETS OF ROIs

%X = sData.stats.sessionAvs.plotXAxis;
binSize = sData.behavior.trialMatrices.meta.binSize;
%nPlaceCells = sData.imdata.placeCells.nPlaceCells;
rewardZone = sData.behavior.rewardZone;
viewDistance = sData.behavior.viewDistance;
smoothSpan = 10;

binnedRois = sData.imdata.binnedRoisDeconv;

hitTrials = nansum(sData.behavior.trialMatrices.rewardInBin(:,1:(120+rewardZone+20)/binSize),2);
hitTrials(hitTrials >= 1) = 1; % Correct if there were multiple reward

missTrials = find(hitTrials==0);
hitTrials = find(hitTrials);



trialSubset1 = intersect(hitTrials,sData.trials.NBCTrials); 
trialSubset2 = intersect(missTrials,sData.trials.NBCTrials);

% correct for different n numbers
trialSubset3 = trialSubset1(randi(numel(trialSubset1),[1 numel(trialSubset2)]));
trialSubset1 = trialSubset3;

% Keap the same, usually BCN light off!
dataForSorting = sData.imdata.avBinnedRois.avBinnedRoisBCNLightOffDeconv;





binNumber = size(binnedRois,2);
nROIs = size(binnedRois,3);

%%%% create position tuning curves
posTuningCurves1(nROIs,binNumber) = NaN;
for roi = 1:1:nROIs 
    posTuningCurves1(roi,:) = nanmean(binnedRois(trialSubset1,:,roi));
end

posTuningCurves2(nROIs,binNumber) = NaN;
for roi = 1:1:nROIs 
    posTuningCurves2(roi,:) = nanmean(binnedRois(trialSubset2,:,roi));
end


%inputData = posTuningCurves2;

subset1 = vr.sortROIs(dataForSorting,sData.imdata.placeCells.placeCells);
%subset2 = vr.sortROIs(dataForSorting,PCsOnlyInstim);
%subset3 = vr.sortROIs(dataForSorting,PCsInBoth);
%subset4 = vr.sortROIs(dataForSorting,notPCsInEither);


data1 = posTuningCurves1(subset1,:);
data2 = posTuningCurves2(subset1,:);
%data3 = inputData(subset3,:);
%data4 = inputData(subset4,:);

data1 = smoothdata(data1,2,'gaussian',smoothSpan);
data2 = smoothdata(data2,2,'gaussian',smoothSpan);
%data3 = smoothdata(data3,2,'gaussian',smoothSpan);
%data4 = smoothdata(data4,2,'gaussian',smoothSpan);




Xax = (-120+binSize):binSize:(binNumber*binSize)-120;
cMax = quantile(data1(data1>0),0.95); %0.02;
cMin = 0;
nPlaceCells = size(data1,1);


figure('Color','white','Position',[0 0 800 800])
%set(0,'DefaultAxesFontSize',12);

subplot(2,2,1);

imagesc(Xax,nPlaceCells:1,data1)
hold on
rectangle('Position',[0,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone-viewDistance,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+20,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+80,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');

%title('\fontsize{13}Light OFF trials')
%t = title('\fontsize{14}Light OFF trials');
%t.FontWeight = 'normal';
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
%xlabel('position in Unity (virtual cm)');
y = ylabel(['\fontsize{12}sorted ROIs']);
%y.FontWeight = 'bold';
c.Label.String = '\fontsize{12}deconvolved activity (a.u.)';
%clabel('DFF (%)');
caxis([cMin cMax]);


subplot(2,2,2);

imagesc(Xax,nPlaceCells:1,data2)
hold on
rectangle('Position',[0,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone-viewDistance,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+20,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+80,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');

%t = title('\fontsize{14}Light ON trials');
%t.FontWeight = 'normal';
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
%xlabel('position in Unity (virtual cm)');
%ylabel('sorted ROIs');
c.Label.String = '\fontsize{12}deconvolved activity (a.u.)';
%clabel('DFF (%)');
caxis([cMin cMax]);


subplot(2,2,3);

imagesc(Xax,nPlaceCells:1,data3)
hold on
rectangle('Position',[0,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone-viewDistance,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+20,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+80,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');

%title('Place cells from all trials DFF')
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('position in Unity (virtual cm)');
y = ylabel(['\fontsize{12}sorted ROIs']);
c.Label.String = '\fontsize{12}deconvolved activity (a.u.)';
%clabel('DFF (%)');
caxis([cMin cMax]);


subplot(2,2,4);

imagesc(Xax,nPlaceCells:1,data4)
hold on
rectangle('Position',[0,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone-viewDistance,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+20,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZone+80,0,2,nPlaceCells],'FaceColor',[1 1 1],'EdgeColor','none');

%title('Place cells from all trials DFF')
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('position in Unity (virtual cm)');
%ylabel('sorted ROIs');
c.Label.String = '\fontsize{12}deconvolved activity (a.u.)';
%clabel('DFF (%)');
caxis([cMin cMax]);


%saveas(gcf,strcat(fullfile(sDataDir,'sortedPlaceCellsDeconv'),'.png'));

%               %               %
%   %       %   %   %       %   %   %
%   %   %   %   %   %   %   %   %   %   %   
%   %       %   %   %       %   %   %
%               %               %

% Plot summary projection

nLightOffTrialsBCN = sData.stats.trialNumbers.nLightOffTrialsBCN;
nLightOnTrialsBCN = sData.stats.trialNumbers.nLightOnTrialsBCN;
nLightOffTrialsNBC = sData.stats.trialNumbers.nLightOffTrialsNBC;
nLightOnTrialsNBC = sData.stats.trialNumbers.nLightOnTrialsNBC;

sumProj(1,:) = normalize(nanmean(data1)); % normalize() handles NaNs in data (not like zscore()) and the defailt method is zscoring
sumProj(2,:) = normalize(nanmean(data2));
%sumProj(3,:) = normalize(nanmean(data3));
%sumProj(4,:) = normalize(nanmean(data4));

%sumStd(1,:) = std(data1);
%sumStd(2,:) = std(data2);
%sumStd(3,:) = std(data3);
%sumStd(4,:) = std(data4);


% Plot summary figure
figure('Color','white')
Xmin = min(Xax);
Xmax = max(Xax);
Ymin = min(sumProj(:))*1.1;
Ymax = max(sumProj(:))*1.3;


hold on
%rectangle('Position',[-1,0,2,nPlaceCells],'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
%rectangle('Position',[rewardZone-viewDistance,0,2,nPlaceCells],'FaceColor',[1 0.8 0.8],'EdgeColor','none');
%rectangle('Position',[rewardZone,0,2,nPlaceCells],'FaceColor',[0.8 1 0.8],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nPlaceCells],'FaceColor',[0.7 1 0.7],'EdgeColor','none');
%rectangle('Position',[rewardZone+79,0,2,nPlaceCells],'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');

rectangle('Position',[-118*0.99,Ymin*0.99,118,-Ymin+Ymax],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,Ymin*0.99,20,-Ymin+Ymax],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),Ymin*0.99,0.5,-Ymin+Ymax],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
rectangle('Position',[rewardZone+80,Ymin*0.99,120,-Ymin+Ymax],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');

plot(Xax,sumProj')
xlabel('\fontsize{12}Position in unity (virtual cm)');
ylabel('\fontsize{12}Mean deconvolved activity (Z-score)');
xlim([Xmin Xmax]);
ylim([Ymin Ymax]);




legend(strcat('Hit trials (n=',num2str(numel(trialSubset1)),')'),...
    strcat('Miss trials (n=',num2str(numel(trialSubset2)),')'),...
   'Location','north')



legend(strcat('PCs only in ctrl (n=',num2str(numel(PCsOnlyInCtrl)),')'),...
    strcat('PCs only in stim. (n=',num2str(numel(PCsOnlyInstim)),')'),...
    strcat('PCs in both (n=',num2str(numel(PCsInBoth)),')'),...
    strcat('Not PCs in either (n=',num2str(numel(notPCsInEither)),')'),'Location','north')





legend(strcat('Beaconed / light-off (n=',num2str(nLightOffTrialsBCN),')'),...
    strcat('Beaconed / light-on (n=',num2str(nLightOnTrialsBCN),')'),...
    strcat('Non-beaconed / light-off (n=',num2str(nLightOffTrialsNBC),')'),...
    strcat('Non-beaconed / light-on (n=',num2str(nLightOnTrialsNBC),')'),'Location','north')

%saveas(gcf,strcat(fullfile(sDataDir,'summaryProjection'),'.png'));

% A(1,:) = min(data1);
% A(2,:) = nanmean(data1);
% A(3,:) = max(data1);

% plotshaded(Xax,A,'r');
%lines(1)



%
%
%                    /\
%                   //\\
%                  ///\\\
%                 ///||\\\
%                /// || \\\
%               //   ||   \\
%                    ||
%                    ||
%                    ||
%

