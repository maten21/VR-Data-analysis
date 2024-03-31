
function fig = plotSortedRoiMatrixEvenOdd(sData)


% Select context to align here:
    t = 1;
    fov = 1;
    subsetOfInterest = sData.imdata.activeROIs;
    
    trials = sData.trials.contextsMeta(t).trials;

    oddTrials = trials(1:2:numel(trials));
    evenTrials = trials(2:2:numel(trials));
    
    
    data = normalize(permute(sData.imdata(fov).binnedRoisDff,[3 2 1]),2);
    % data = permute(sData.imdata(fov).binnedRoisDeconvRate,[3 2 1]);
    
    smoothSpan = 5;    
    posTuningCurves = nanmean(data(:,:,evenTrials),3);
    sortedROIs = vr.sortROIs(posTuningCurves,subsetOfInterest,smoothSpan);
   % sortedROIs = vr.sortROIs(posTuningCurves); %,subsetOfInterest,smoothSpan);
    nSortedROIs = numel(sortedROIs);
    
 

% Define matrices for subplots 

    %trials = oddTrials;

M1 = nanmean(data(sortedROIs,:,evenTrials),3);
M2 = nanmean(data(sortedROIs,:,oddTrials),3);
Xax = sData.stats.sessionAvs.plotXAxis;
cMin = -0.25;
cMax = 2;

fig = figure('Color','white','Position',[0 0 800 400]);
%set(0,'DefaultAxesFontSize',12);
cLabel = 'DFF (Z-score)';
% cLabel = 'Deconv. activity rate (Z-score)';

subplot(1,2,1)

imagesc(Xax,nSortedROIs:1,M1)
hold on
%rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

% rectangle('Position',[rewardZones(1,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(1,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');

%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

%title('\fontsize{13}Light OFF trials')
t = title('Even trials');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([cMin cMax]);


subplot(1,2,2)

imagesc(Xax,nSortedROIs:1,M2)
hold on
%rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

% rectangle('Position',[rewardZones(1,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(1,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');

%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

%title('\fontsize{13}Light OFF trials')
t = title('Odd trials');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([cMin cMax]);
   





fig = figure('Color','white','Position',[0 0 500 400]);
%set(0,'DefaultAxesFontSize',12);
cLabel = 'Corr. coef.';
% cLabel = 'Deconv. activity rate (Z-score)';

%subplot(1,2,1)

imagesc(Xax,Xax,corr(M1,M2))
hold on
%rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

% rectangle('Position',[rewardZones(1,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(1,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');

%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

%title('\fontsize{13}Light OFF trials')
t = title('Even trials');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([-0.5 1]);


% figure

hold on
plot(Xax,max(corr(M1,M2)))




end