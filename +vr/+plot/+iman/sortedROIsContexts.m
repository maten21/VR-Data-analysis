%% PLOT sorted ROI matrices with separated trial types DFF
%X = sData.stats.sessionAvs.plotXAxis;

% sort ROIs

    % posTunCurves = sData.imdata.avBinnedRois.avBinnedRoisDeconv{1};    
    % sortedROIs = vr.sortROIs(posTunCurves,subsetOfInterest,smoothSpan);    
    
    % Down sample control trials 
%    n = sData.trials.trialTypesMeta(1).nTrials;
%    k = sData.trials.trialTypesMeta(2).nTrials; % to randomly select the same number of trials   
%    ctrlTrials = sData.trials.trialTypesMeta(1).trials(randperm(n,k)); %returns a row vector containing k unique integers selected randomly from 1 to n inclusive
%    complementerTrials = setdiff(sData.trials.trialTypesMeta(1).trials,ctrlTrials);    
    
%%

sData0 = sData;

sData.imdata = sData0.imdata(1);
%sData.trials = sData0.trials;




%% PLOT Trial Types DFF
%    threshold = 0.3;
%    freqThreshold = 1.5;
    
    % subsetOfInterest = find(activeTrialFraction > threshold);
    % subsetOfInterest = intersect(find(activeTrialFractionCtrl > threshold),find(peakFreqs > freqThreshold));
    % numel(subsetOfInterest)
    nROIs = sData.imdata.nROIs;
    
     subsetOfInterest = 1:1:nROIs;
    
    % active ROIs in any of the trial types
%    subsetOfInterest = intersect(union(find(activeTrialFractionInBlock > threshold),union(find(activeTrialFractionHB > threshold),find(activeTrialFractionOTP > threshold))),...
%        union(find(peakFreqsCtrl > freqThreshold),union(find(peakFreqsHB > freqThreshold),find(peakFreqsOTP > freqThreshold))));

activity = [sData.imdata.roiStat.activityLevel];
subsetOfInterest = sData.imdata.activeROIs;

% Place cell in either block:
subsetOfInterest = union(union(sData.imdata.placeCells(1).placeCells,sData.imdata.placeCells(2).placeCells),...
    union(sData.imdata.placeCells(3).placeCells,sData.imdata.placeCells(4).placeCells));


%subsetOfInterest = sData.imdata.inactiveROIs;
%subsetOfInterest = sData.imdata.tunedROIs;
%subsetOfInterest = sData.imdata.untunedROIs;

%subsetOfInterest = intersect(inactiveROIs,untunedROIs);    
%subsetOfInterest = intersect(inactiveROIs,tunedROIs);
%subsetOfInterest = intersect(activeROIs,tunedROIs);
%subsetOfInterest = intersect(activeROIs,untunedROIs);



    % subsetOfInterest = union(find(activeTrialFractionCtrl > threshold),union(find(activeTrialFractionHB > threshold),find(activeTrialFractionOTP > threshold)));
    % numel(subsetOfInterest)
    
    % subsetOfInterestComplementer = setdiff(1:1:nROIs,subsetOfInterest);
    % subsetOfInterest = subsetOfInterestComplementer;
    % subsetOfInterest = sData.imdata.placeCells.placeCells;  
    
    
    
%nFOVs = length(sData.imdata);    
    
    
    
    trials = sData.trials.contextsMeta(1).trials;
    if rem(numel(trials),2) > 0
        trials = [trials trials(end)];
    end
    
    A = reshape(trials,[2 ceil(numel(trials)/2)]);
    oddTrials = A(1,:);
    evenTrials = A(2,:);
    
    
    data = normalize(permute(sData.imdata.binnedRoisDff,[3 2 1]),2);
    % data = permute(sData.imdata.binnedRoisDeconvRate,[3 2 1]);
    
    smoothSpan = 5;    
    posTuningCurves = nanmean(data(:,:,evenTrials),3);
    sortedROIs = vr.sortROIs(posTuningCurves,subsetOfInterest,smoothSpan);
   % sortedROIs = vr.sortROIs(posTuningCurves); %,subsetOfInterest,smoothSpan);
    nSortedROIs = numel(sortedROIs);
    
 

% Define matrices for subplots 
trials = sData.trials.contextsMeta(1).trials;
% trials = oddTrials;
M1 = nanmean(data(sortedROIs,:,trials),3);

trials = sData.trials.contextsMeta(2).trials;
% trials = oddTrials;
M2 = nanmean(data(sortedROIs,:,trials),3);

trials = sData.trials.contextsMeta(3).trials;
M3 = nanmean(data(sortedROIs,:,trials),3);

trials = sData.trials.contextsMeta(4).trials;
M4 = nanmean(data(sortedROIs,:,trials),3);

    
    
%{    
    data = permute(sData.imdata.binnedRoisDff,[3 2 1]);
    % data = permute(sData.imdata.binnedRoisDeconvRate,[3 2 1]);
    
    smoothSpan = 5;    
    sortedROIs = vr.sortROIs(nanmean(data(:,:,complementerTrials),3),subsetOfInterest,smoothSpan);
    nSortedROIs = numel(sortedROIs);
    
 

% Define matrices for subplots 
trials = sData.trials.contextsMeta(1).trials;
% trials = ctrlTrials;
M1 = normalize(nanmean(data(sortedROIs,:,trials),3),2);

trials = sData.trials.contextsMeta(2).trials;
% trials = ctrlTrials;
M2 = normalize(nanmean(data(sortedROIs,:,trials),3),2);

trials = sData.trials.contextsMeta(3).trials;
M3 = normalize(nanmean(data(sortedROIs,:,trials),3),2);

trials = sData.trials.contextsMeta(4).trials;
M4 = normalize(nanmean(data(sortedROIs,:,trials),3),2);
%}

binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;

rew = mean(sData.behavior.trialMatrices.rewardInBin);
rew(isnan(rew)) = 0;
rewardZones = sData.stats.sessionAvs(1).plotXAxis(find(rew));
rewardZones = rewardZones - binSize;

corridorLength = binNumber*binSize;
viewDistance = 50;

Xax = sData.stats.sessionAvs(1).plotXAxis;
cMin = -0.25;
cMax = 2;
% nSortedROIs = numel(sortedROIs);


figure('Color','white','Position',[0 0 800 800])
%set(0,'DefaultAxesFontSize',12);
cLabel = 'DFF (Z-score)';
% cLabel = 'Deconv. activity rate (Z-score)';

subplot(2,2,1);

imagesc(Xax,nSortedROIs:1,M1)
hold on
rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(3),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

%title('\fontsize{13}Light OFF trials')
t = title('Block 1 - Context A');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
%xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([cMin cMax]);


subplot(2,2,2);

imagesc(Xax,nSortedROIs:1,M2)
hold on
rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(3),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

t = title('Block 2 - Context B');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
%xlabel('\fontsize{13}Track position (cm)');
%ylabel('sorted ROIs');
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([cMin cMax]);


subplot(2,2,3);

imagesc(Xax,nSortedROIs:1,M3)
hold on
rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(3),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

t = title('Block 3 - Context A');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([cMin cMax]);


subplot(2,2,4);

imagesc(Xax,nSortedROIs:1,M4)
hold on
rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(3),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

t = title('Block 4 - Context B');
t.FontSize = 14;
t.FontWeight = 'normal';
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
%y = ylabel('\fontsize{12}sorted ROIs');
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([cMin cMax]);

suptitle('All are sorted based on Block 1');

saveas(gcf,strcat(fullfile(sDataDir,'sortedROIsDff_1_allCells_FOV4'),'.png'));

