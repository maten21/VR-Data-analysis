%% PLOT sorted ROI matrices with separated trial types DFF
%X = sData.stats.sessionAvs.plotXAxis;

% sort ROIs

    % posTunCurves = sData.imdata(fov).avBinnedRois.avBinnedRoisDeconv{1};    
    % sortedROIs = vr.sortROIs(posTunCurves,subsetOfInterest,smoothSpan);    
    
    % Down sample control trials 
%    n = sData.trials.trialTypesMeta(1).nTrials;
%    k = sData.trials.trialTypesMeta(2).nTrials; % to randomly select the same number of trials   
%    ctrlTrials = sData.trials.trialTypesMeta(1).trials(randperm(n,k)); %returns a row vector containing k unique integers selected randomly from 1 to n inclusive
%    complementerTrials = setdiff(sData.trials.trialTypesMeta(1).trials,ctrlTrials);    
    
function sData = remapPlotsRoiMatricesOpto(sData,sDataDir) 


for fov = 1:1:length(sData.imdata)

%% PLOT Trial Types DFF

% subsetOfInterest = 1:1:nROIs;
    
    % active ROIs in any of the trial types
%    subsetOfInterest = intersect(union(find(activeTrialFractionInBlock > threshold),union(find(activeTrialFractionHB > threshold),find(activeTrialFractionOTP > threshold))),...
%        union(find(peakFreqsCtrl > freqThreshold),union(find(peakFreqsHB > freqThreshold),find(peakFreqsOTP > freqThreshold))));

subsetOfInterest = sData.imdata(fov).activeROIs;

% Place cell in either block:
%subsetOfInterest = union(union(sData.imdata(fov).placeCells(1).placeCells,sData.imdata(fov).placeCells(2).placeCells),...
%    union(sData.imdata(fov).placeCells(3).placeCells,sData.imdata(fov).placeCells(4).placeCells));


%subsetOfInterest = sData.imdata(fov).inactiveROIs;
%subsetOfInterest = sData.imdata(fov).tunedROIs;
%subsetOfInterest = sData.imdata(fov).untunedROIs;

%subsetOfInterest = intersect(inactiveROIs,untunedROIs);    
%subsetOfInterest = intersect(inactiveROIs,tunedROIs);
%subsetOfInterest = intersect(activeROIs,tunedROIs);
%subsetOfInterest = intersect(activeROIs,untunedROIs);



binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
rewardZones = [];
optStim = struct; % for rectangle plot

for context = 1:1:length(sData.trials.contextsMeta)
rew = mean(sData.behavior.trialMatrices.rewardInBin(sData.trials.contextsMeta(context).trials,:));
rew(isnan(rew)) = 0;
rewardZ = sData.stats.sessionAvs(1).plotXAxis(find(rew));

rewardZones(context,1) = rewardZ(1); 
rewardZones(context,2) = rewardZ(find(diff(discretize(rewardZ,2)))+1); 

optStim(context).isStim = numel(sData.trials.stimProtocols(context).from) > 0;
optStim(context).from = sData.trials.stimProtocols(context).from;
optStim(context).length = sData.trials.stimProtocols(context).to - sData.trials.stimProtocols(context).from;

end
rewardZones = rewardZones - binSize;





for sortMat = [1, 2, 3, 4, 5]
    

    
    
    % Select context to align here:
    trials = sData.trials.contextsMeta(sortMat).trials;

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
if sortMat == 1
    trials = oddTrials;
else
    trials = sData.trials.contextsMeta(1).trials;
end
M1 = nanmean(data(sortedROIs,:,trials),3);

if sortMat == 2
    trials = oddTrials;
else
    trials = sData.trials.contextsMeta(2).trials;
end
M2 = nanmean(data(sortedROIs,:,trials),3);

if sortMat == 3
    trials = oddTrials;
else
    trials = sData.trials.contextsMeta(3).trials;
end
M3 = nanmean(data(sortedROIs,:,trials),3);

if sortMat == 4
    trials = oddTrials;
else
    trials = sData.trials.contextsMeta(4).trials;
end
M4 = nanmean(data(sortedROIs,:,trials),3);

if sortMat == 5
    trials = oddTrials;
else
    trials = sData.trials.contextsMeta(5).trials;
end
M5 = nanmean(data(sortedROIs,:,trials),3);

    
%{    
    data = permute(sData.imdata(fov).binnedRoisDff,[3 2 1]);
    % data = permute(sData.imdata(fov).binnedRoisDeconvRate,[3 2 1]);
    
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

corridorLength = binNumber*binSize;
viewDistance = 50;

Xax = sData.stats.sessionAvs(1).plotXAxis;
cMin = -0.25;
cMax = 2;
% nSortedROIs = numel(sortedROIs);


figure('Color','white','Position',[0 0 1200 800])
%set(0,'DefaultAxesFontSize',12);
cLabel = 'DFF (Z-score)';
% cLabel = 'Deconv. activity rate (Z-score)';

subplot(2,3,2);

imagesc(Xax,nSortedROIs:1,smoothdata(M1,2,"gaussian",smoothSpan,"omitnan"))
hold on
% rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(1,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(1,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

%title('\fontsize{13}Light OFF trials')
t = title('Block 1 - Fam - No-stim.');
if sortMat ~= 1
t.FontWeight = 'normal';
end
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

% %{
subplot(2,3,3);

imagesc(Xax,nSortedROIs:1,smoothdata(M5,2,"gaussian",smoothSpan,"omitnan"))
hold on
% rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(5,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(5,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

t = title('Block 5 - Fam - Opto','Color','r');
if sortMat ~= 5
t.FontWeight = 'normal';
end
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
% %}
subplot(2,3,4);

imagesc(Xax,nSortedROIs:1,smoothdata(M2,2,"gaussian",smoothSpan,"omitnan"))
hold on
% rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(2,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(2,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

t = title('Block 2 - New - Opto','Color','r');
if sortMat ~= 2
t.FontWeight = 'normal';
end
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


subplot(2,3,5);

imagesc(Xax,nSortedROIs:1,smoothdata(M3,2,"gaussian",smoothSpan,"omitnan"))
hold on
% rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(3,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(3,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

t = title('Block 3 - New - No-stim');
if sortMat ~= 3
t.FontWeight = 'normal';
end
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


subplot(2,3,6);

imagesc(Xax,nSortedROIs:1,smoothdata(M4,2,"gaussian",smoothSpan,"omitnan"))
hold on
% rectangle('Position',[50,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(4,1),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(4,2),0,2,nSortedROIs],'FaceColor',[0 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone-viewDistance,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZone,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+20,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');
%rectangle('Position',[rewardZone+80,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');

t = title('Block 4 - New - Opto','Color','r');
t.FontSize = 14;
if sortMat ~= 4
t.FontWeight = 'normal';
end
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
%y = ylabel('\fontsize{12}sorted ROIs');
c.Label.String = cLabel;
%clabel('DFF (%)');
caxis([cMin cMax]);


suptitle([sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation ' - Block ' num2str(sortMat) ' sorted']);

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation '-sortedROIsDff_' num2str(sortMat) '_activeCells']),'.png')); % 'sortedROIsDff_2_allPlaceCells'

% saveas(gcf,strcat(fullfile(sDataDir,'sortedROIsDffComplementer'),'.png'));
% saveas(gcf,strcat(fullfile(sDataDir,'sortedROIsDffPlaceCells'),'.png'));

% Auto-correlate

% figure;
% imagesc(Xax,Xax,corr(M1))


end

%% Cross corr of trial blocks

% Define matrices for subplots 
trials = sData.trials.contextsMeta(1).trials;
B1 = nanmean(data(subsetOfInterest,:,trials),3);

trials = sData.trials.contextsMeta(2).trials;
B2 = nanmean(data(subsetOfInterest,:,trials),3);

trials = sData.trials.contextsMeta(3).trials;
B3 = nanmean(data(subsetOfInterest,:,trials),3);

trials = sData.trials.contextsMeta(4).trials;
B4 = nanmean(data(subsetOfInterest,:,trials),3);

trials = sData.trials.contextsMeta(5).trials;
B5 = nanmean(data(subsetOfInterest,:,trials),3);




B13 = corr(B1,B3); 
B15 = corr(B1,B5);
B32 = corr(B3,B2);
B34 = corr(B3,B4); 

sData.imdata.correlations.xBlocks.corrCoef.B1xB3 = B13;
sData.imdata.correlations.xBlocks.corrCoef.B1xB5 = B15;
sData.imdata.correlations.xBlocks.corrCoef.B3xB2 = B32;
sData.imdata.correlations.xBlocks.corrCoef.B3xB4 = B34;

sData.imdata.correlations.xBlocks.maxCorrCoef.B1xB3 = max(B13);
sData.imdata.correlations.xBlocks.maxCorrCoef.B1xB5 = max(B15);
sData.imdata.correlations.xBlocks.maxCorrCoef.B3xB2 = max(B32);
sData.imdata.correlations.xBlocks.maxCorrCoef.B3xB4 = max(B34);

lightBlue = [0.3010 0.7450 0.9330];

figure('Color','white','Position',[0 0 1200 800])
hold on


subplot(2,3,3)
hold on

rectangle('Position',[optStim(2).from+2, 0.96, optStim(2).length, 0.03],'FaceColor',[1 0.2 0.2],'EdgeColor','none');

rectangle('Position',[200, 0.9, 25, 0.1],'FaceColor','k','EdgeColor','none');
rectangle('Position',[225, 0.9, 25, 0.1],'FaceColor','g','EdgeColor','none');
rectangle('Position',[200, 0.8, 25, 0.1],'FaceColor',lightBlue,'EdgeColor','none');
rectangle('Position',[225, 0.8, 25, 0.1],'FaceColor','b','EdgeColor','none');

plot(Xax,smoothdata(max(B13),'gaussian',smoothSpan),'k-','LineWidth',2)
plot(Xax,smoothdata(max(B15),'gaussian',smoothSpan),'g-','LineWidth',2)
plot(Xax,smoothdata(max(B32),'gaussian',smoothSpan),'Color',lightBlue,'LineWidth',2)
plot(Xax,smoothdata(max(B34),'gaussian',smoothSpan),'b-','LineWidth',2)

xlabel('\fontsize{13}Track position (cm)');
ylabel('\fontsize{13}Max corr. coef.'); 

ylim([0 1])
% tit = title('Max corr. coef.');
%tit.FontWeight = 'normal';
%tit.FontSize = 14;
% tit.Color = 'r';
% legend(subPlotTitles,'Location','southeast')




subplot(2,3,1);
imagesc(Xax,Xax,B13)

colormap(gca,jet);
caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
% t = title('Block 1 - Context A'); t.FontSize = 15; t.FontWeight = 'normal';
ylabel(['\fontsize{16}Block 1 - Fam. - No-stim.' newline '\fontsize{13}Track position (cm)'])
%xlabel(['\fontsize{16}Block 3 - New - No-stim.' newline '\fontsize{13}Track position (cm)'])
xlabel([ '\fontsize{13}Track position (cm)'])
tit = title('\fontsize{16}Block 3 - New - No-stim.');
tit.FontWeight = 'normal';

subplot(2,3,2);
hold on

imagesc(Xax,Xax,B15)

rectangle('Position',[optStim(2).from,optStim(2).from,optStim(2).length,optStim(2).length],'FaceColor','none','EdgeColor',[1 0.2 0.2],'LineWidth',3);

colormap(gca,jet);
caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
% t = title('Block 2 - Context B'); t.FontSize = 15; t.FontWeight = 'normal';
%ylabel(['\fontsize{16}Block 1 - Fam. - No-stim.' newline '\fontsize{13}Track position (cm)'])
%xlabel(['\fontsize{16}Block 5 - Fam.- Opto.' newline '\fontsize{13}Track position (cm)'])
xlabel(['\fontsize{13}Track position (cm)'])

tit = title('\fontsize{16}Block 5 - Fam.- Opto.');
tit.FontWeight = 'normal';
tit.Color = 'r';

subplot(2,3,4);
hold on

imagesc(Xax,Xax,B32)

rectangle('Position',[optStim(2).from,optStim(2).from,optStim(2).length,optStim(2).length],'FaceColor','none','EdgeColor',[1 0.2 0.2],'LineWidth',3);

colormap(gca,jet);
caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
ylabel(['\fontsize{16}Block 3 - New - No-stim.' newline '\fontsize{13}Track position (cm)'])
%xlabel(['\fontsize{16}Block 2 - New - Opto.' newline '\fontsize{13}Track position (cm)'])
xlabel(['\fontsize{13}Track position (cm)'])

tit = title('\fontsize{16}Block 2 - New first time - Opto.');
tit.FontWeight = 'normal';
tit.Color = 'r';


subplot(2,3,5);
hold on
imagesc(Xax,Xax,B34)

rectangle('Position',[optStim(2).from,optStim(2).from,optStim(2).length,optStim(2).length],'FaceColor','none','EdgeColor',[1 0.2 0.2],'LineWidth',3);

colormap(gca,jet);
caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
%ylabel(['\fontsize{16}Block 3 - New - No-stim.' newline '\fontsize{13}Track position (cm)'])
%xlabel(['\fontsize{16}Block 4 - New - Opto.' newline '\fontsize{13}Track position (cm)'])
xlabel(['\fontsize{13}Track position (cm)'])

tit = title('\fontsize{16}Block 4 - New - Opto.');
tit.FontWeight = 'normal';
tit.Color = 'r';



suptitle([sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation ]) %newline ...
%    'Population vector cross correlation of trial blocks']);

% suptitle('\fontsize{16}Population vector cross correlation of trial blocks');

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation '-crossCorrDff']),'.png'));



%% Diagonal of cross corrs
%{
shift = floor(numel(Xax)/2);
tempM = B13;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
Cs1 = tempM;

tempM = B15;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
Cs2 = tempM;

tempM = B32;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
Cs3 = tempM;

tempM = B34;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
Cs4 = tempM;
%{
tempM = C5;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
Cs5 = tempM;

tempM = C6;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
Cs6 = tempM;
%}



figure('Color','white','Position',[0 0 800 330])
suptitle(['\fontsize{13}' sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation newline]);

subplot(1,2,1)
hold on
% rectangle('Position',[50,-0.39,140,1.39],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
%rectangle('Position',[49,-0.4,2,1.4],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
%rectangle('Position',[189,-0.4,2,1.4],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
rectangle('Position',[rewardZones(1,1)-2,0.95,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(1,2)-2,0.95,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(2,1)-2,0.9,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
rectangle('Position',[rewardZones(2,2)-2,0.9,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
plot(Xax,[Cs1(:,shift) Cs2(:,shift) Cs3(:,shift) Cs4(:,shift) ]) % Cs5(:,shift) Cs6(:,shift)])

xlabel('\fontsize{13}Track position (cm)')
ylabel('\fontsize{13}Corr. coefficient at diagonal')
%legend('A1 X A2','B1 X A2','A1 X B2','B1 X B2')
text(110,1.05,'Reward positions','color',[0 1 1])
% text(2,-0.3,'Identical')
% text(200,-0.3,'Identical')
% text(70,-0.3,'Unique VR context')
ax = gca;
ax.TickDir = 'out';
ylim([-0.4 1])


%{
subplot(1,2,2)
hold on
bins = [1:24 96:125];
plot((-shift:shift)*2,nanmean(Cs1(bins,:)))
plot((-shift:shift)*2,nanmean(Cs4(bins,:)))
plot((-shift:shift)*2,nanmean(Cs2(bins,:)))
plot((-shift:shift)*2,nanmean(Cs3(bins,:)))
plot((-shift:shift)*2,nanmean(Cs5(bins,:)))
plot((-shift:shift)*2,nanmean(Cs6(bins,:)))
text(-40,1.05,'Identical part','color',[0 0 0])
%{
subplot(1,2,2)
hold on
plot((-shift:shift)*2,nanmean(Cs1))
plot((-shift:shift)*2,nanmean(Cs2))
plot((-shift:shift)*2,nanmean(Cs3))
plot((-shift:shift)*2,nanmean(Cs4))
%}


xlabel('\fontsize{13}Distance from diagonal (cm)')
ylabel('\fontsize{13}Mean corr. coefficient')
legend('A1 X A2','B1 X B2','B1 X A2','A1 X B2','A1 X B1','A2 X B2')
ax = gca;
ax.TickDir = 'out';
%}

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation '-crossCorrCoefDiagonalDff']),'.png'));


%{

figure('Color','white','Position',[0 0 400 300])
hold on

mymap = lines;
bins = [1:24 96:125];
plotMeanStd((-shift:shift)*2, [Cs2(bins,:)', Cs3(bins,:)', Cs5(bins,:)', Cs6(bins,:)'],mymap(1,:),1)
bins = 25:95;
plotMeanStd((-shift:shift)*2, [Cs2(bins,:)', Cs3(bins,:)', Cs5(bins,:)', Cs6(bins,:)'],mymap(2,:),1)

xlabel('\fontsize{13}Distance from diagonal (cm)')
ylabel('\fontsize{13}Mean corr. coefficient')
% legend('Identical visual context','Changing visual context')
title( [sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation newline 'Identical vs changing visual context'])

saveas(gcf,strcat(fullfile(sDataDir,['Fov' num2str(fov) '-' sData.imdata(fov).fovLocation '-crossCorrIdenticalVsUniqueDff']),'.png'));

%}


%}

%% correlate sorting position

%{
smoothSpan = 5;  

      
[~, sortIndexes1] = vr.sortROIs(nanmean(data(:,:,sData.trials.contextsMeta(1).trials),3),subsetOfInterest,smoothSpan);
[~, sortIndexes2] = vr.sortROIs(nanmean(data(:,:,sData.trials.contextsMeta(2).trials),3),subsetOfInterest,smoothSpan);

[fitX, fitY, R2, P] = linFit(sortIndexes1,sortIndexes2); %[fitX, fitY, R2, P, a, b] = linFit(sortIndexes1,sortIndexes2);

figure
hold on
plot(sortIndexes1,sortIndexes2,'.')
plot(fitX, fitY,'-')
%}


%% Auto correlation Odd X Even



figure('Color','white','Position',[0 0 1200 800]);
subPlotPos = [2, 4, 5, 6, 3];
subPlotTitles = {'Block 1 - Fam - No-stim.', 'Block 2 - New - Opto.', 'Block 3 - New - No-stim.', 'Block 4 - New - Opto.', 'Block 5 - Fam - Opto.'};

myMap = lines;


    fov = 1;
    subsetOfInterest = sData.imdata.activeROIs;
    
    
    data = normalize(permute(sData.imdata(fov).binnedRoisDff,[3 2 1]),2);
    % data = permute(sData.imdata(fov).binnedRoisDeconvRate,[3 2 1]);
    
    smoothSpan = 5;    

for t = 1:1:5
% Select context to align here:
    %t = 5;

    trials = sData.trials.contextsMeta(t).trials;
    oddTrials = trials(1:2:numel(trials));
    evenTrials = trials(2:2:numel(trials));
    
    posTuningCurves = nanmean(data(:,:,evenTrials),3);
    sortedROIs = vr.sortROIs(posTuningCurves,subsetOfInterest,smoothSpan);
   % sortedROIs = vr.sortROIs(posTuningCurves); %,subsetOfInterest,smoothSpan);
    nSortedROIs = numel(sortedROIs);
    
 

% Define matrices for subplots 

    %trials = oddTrials;

M1 = nanmean(data(sortedROIs,:,evenTrials),3);
M2 = nanmean(data(sortedROIs,:,oddTrials),3);
Xax = sData.stats.sessionAvs.plotXAxis;

sData.imdata.correlations.xEvenOdd(t).maxCorrCoef = max(corr(M1,M2));
sData.imdata.correlations.xEvenOdd(t).corrCoef = corr(M1,M2);

subplot(2,3,1)
hold on
rectangle('Position',[optStim(2).from+2, 0.96, optStim(2).length, 0.03],'FaceColor',[1 0.2 0.2],'EdgeColor','none');

rectangle('Position',[200, 0.2, 25, 0.1],'FaceColor',myMap(1,:),'EdgeColor','none');
rectangle('Position',[175, 0.1, 25, 0.1],'FaceColor',myMap(2,:),'EdgeColor','none');
rectangle('Position',[200, 0.1, 25, 0.1],'FaceColor',myMap(3,:),'EdgeColor','none');
rectangle('Position',[225, 0.1, 25, 0.1],'FaceColor',myMap(4,:),'EdgeColor','none');
rectangle('Position',[225, 0.2, 25, 0.1],'FaceColor',myMap(5,:),'EdgeColor','none');

plot(Xax,smoothdata(max(corr(M1,M2)),'gaussian',smoothSpan),'LineWidth',2)

xlabel('\fontsize{13}Track position (cm)');
ylabel('\fontsize{13}Max corr. coef.'); 

ylim([0 1])
% tit = title('Max corr. coef.');
%tit.FontWeight = 'normal';
%tit.FontSize = 14;
% tit.Color = 'r';
% legend(subPlotTitles,'Location','southeast')




subplot(2,3,subPlotPos(t))

imagesc(Xax,Xax,corr(M1,M2))
hold on
if optStim(t).isStim
rectangle('Position',[optStim(t).from,optStim(t).from,optStim(t).length,optStim(t).length],'FaceColor','none','EdgeColor',[1 0.2 0.2],'LineWidth',3);
end

ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';

%rectangle('Position',[190,0,2,nSortedROIs],'FaceColor',[1 1 1],'EdgeColor','none');


if optStim(t).isStim
    tit = title(subPlotTitles{t},'Color','r');
else
    tit = title(subPlotTitles{t});
end

tit.FontWeight = 'normal';
tit.FontSize = 14;

colormap(gca,jet);
%c = colorbar;
%c.Label.String = cLabel;
%caxis([-0.5 1]);
%cLabel = 'Corr. coef.';

ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
if subPlotPos(t) == 4
    y = ylabel('\fontsize{13}Track position (cm)'); 
end

end
suptitle(sData.sessionInfo.sessionID)

saveas(gcf,strcat(fullfile(sDataDir,['Fov' num2str(fov) '-' sData.imdata(fov).fovLocation '-autoCorrTrialBlocksOddEvenDff']),'.png'));



vr.plot.iman.mapStabilityQC(sData,sDataDir);







end


end