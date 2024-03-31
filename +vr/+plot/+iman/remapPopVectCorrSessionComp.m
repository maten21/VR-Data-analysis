% Session Comparison population vector correlations

clear
[sDataFiles, filePath] = vr.loadData('light');


% remapPopVectCorr 




% create FOV library
fovLib = struct();
i = 1;
for f = 1:1:length(sDataFiles)

for fov = 1:1:length(sDataFiles{f}.imdata)

    fovLib(i).fileInd = f;
    fovLib(i).fovInd = fov;
    fovLib(i).fovLoc = sDataFiles{f}.imdata(fov).fovLocation;
    
    if strcmp(sDataFiles{f}.imdata(fov).fovLocation(1:2),'RS')
        fovLib(i).isRSC = true;
    else
        fovLib(i).isRSC = false;
    end
    if strcmp(sDataFiles{f}.imdata(fov).fovLocation(1:2),'M2')
        fovLib(i).isM2 = true;
    else
        fovLib(i).isM2 = false;
    end
    if strcmp(sDataFiles{f}.imdata(fov).fovLocation(1:2),'V2')
        fovLib(i).isV2 = true;
    else
        fovLib(i).isV2 = false;
    end
    if strcmp(sDataFiles{f}.imdata(fov).fovLocation(1:2),'PP')
        fovLib(i).isPPC = true;
    else
        fovLib(i).isPPC = false;
    end
    if strcmp(sDataFiles{f}.imdata(fov).fovLocation(1:2),'HP')
        fovLib(i).isHPC = true;
    else
        fovLib(i).isHPC = false;
    end
    
   i = i + 1;
end


end


shiftA1A2 = [];
shiftB1B2 = [];
shiftAB = [];
velA1 = [];
velA2 = [];
velB1 = [];
velB2 = [];
velA = [];
velB = [];
lickA1 = [];
lickA2 = [];
lickB1 = [];
lickB2 = [];
lickA = [];
lickB = [];

for f = 1:1:length(fovLib)

    fInd = fovLib(f).fileInd;
    fov = fovLib(f).fovInd;
    
    roiSubset = sDataFiles{fInd}.imdata(fov).activeROIs;
    
    %data = normalize(permute(sDataFiles{fInd}.imdata(fov).binnedRoisDff,[3 2 1]),2); % Does it need to be normalized for correlation? DFF or Deconv?
     data = permute(sDataFiles{fInd}.imdata(fov).binnedRoisDff,[3 2 1]); 
    %data = permute(sDataFiles{fInd}.imdata(fov).binnedRoisDeconvRate,[3 2 1]); 
    velData = sDataFiles{fInd}.behavior.trialMatrices.binVel;
    lickData = sDataFiles{fInd}.behavior.trialMatrices.lickFreqInBin;
    
% Define matrices for subplots 
trials = sDataFiles{fInd}.trials.contextsMeta(1).trials;
A1 = nanmean(data(roiSubset,:,trials),3);
velA1(:,f) = nanmean(velData(trials,:));
lickA1(:,f) = nanmean(lickData(trials,:));

trials = sDataFiles{fInd}.trials.contextsMeta(2).trials;
B1 = nanmean(data(roiSubset,:,trials),3);
velB1(:,f) = nanmean(velData(trials,:));
lickB1(:,f) = nanmean(lickData(trials,:));

trials = sDataFiles{fInd}.trials.contextsMeta(3).trials;
A2 = nanmean(data(roiSubset,:,trials),3);
velA2(:,f) = nanmean(velData(trials,:));
lickA2(:,f) = nanmean(lickData(trials,:));

trials = sDataFiles{fInd}.trials.contextsMeta(4).trials;
B2 = nanmean(data(roiSubset,:,trials),3);
velB2(:,f) = nanmean(velData(trials,:));
lickB2(:,f) = nanmean(lickData(trials,:));

trials = union(sDataFiles{fInd}.trials.contextsMeta(1).trials,sDataFiles{fInd}.trials.contextsMeta(3).trials);
A = nanmean(data(roiSubset,:,trials),3);
velA(:,f) = nanmean(velData(trials,:));
lickA(:,f) = nanmean(lickData(trials,:));

trials = union(sDataFiles{fInd}.trials.contextsMeta(2).trials,sDataFiles{fInd}.trials.contextsMeta(4).trials);
B = nanmean(data(roiSubset,:,trials),3);
velB(:,f) = nanmean(velData(trials,:));
lickB(:,f) = nanmean(lickData(trials,:));


A1A2 = corr(A1,A2); 
B1B2 = corr(B1,B2);
AB = corr(A,B);




Xax = sDataFiles{fInd}.stats.sessionAvs(1).plotXAxis; 
shift = floor(numel(Xax)/2);
tempM = A1A2;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
shiftA1A2(:,:,f) = tempM;

tempM = B1B2;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
shiftB1B2(:,:,f) = tempM;

tempM = AB;
for i = 1:1:numel(Xax)
    temp = [tempM(i,i:numel(Xax)) tempM(i,1:i-1)];
    tempM(i,:) = [temp((shift+1):numel(Xax)) temp(1:shift)];
end
shiftAB(:,:,f) = tempM;


end


%{

figure('Color','white','Position',[0 0 800 800])
hold on

subplot(2,2,1);
imagesc(Xax,Xax,C1)

% identical context borders:
rectangle('Position',[49.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,49.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

rectangle('Position',[189.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,189.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
% t = title('Block 1 - Context A'); t.FontSize = 15; t.FontWeight = 'normal';
ylabel(['\fontsize{16}Block 3 - Context A' newline '\fontsize{13}Track position (cm)'])
% xlabel('\fontsize{13}Track position (cm)')

subplot(2,2,2);
imagesc(Xax,Xax,C2)

% identical context borders:
rectangle('Position',[49.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,49.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

rectangle('Position',[189.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,189.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
% t = title('Block 2 - Context B'); t.FontSize = 15; t.FontWeight = 'normal';
% xlabel('\fontsize{13}Track position (cm)')

subplot(2,2,3);
imagesc(Xax,Xax,C3)

% identical context borders:
rectangle('Position',[49.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,49.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

rectangle('Position',[189.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,189.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
ylabel(['\fontsize{16}Block 4 - Context B' newline '\fontsize{13}Track position (cm)'])
xlabel(['\fontsize{13}Track position (cm)' newline '\fontsize{16}Block 1 - Context A'])

subplot(2,2,4);
imagesc(Xax,Xax,C4)

% identical context borders:
rectangle('Position',[49.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,49.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

rectangle('Position',[189.5,min(Xax)*0.99,1,-min(Xax)+max(Xax)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[min(Xax)*0.99,189.5,-min(Xax)+max(Xax),1],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box

caxis([-0.3 1])
axis([min(Xax)*1.015 max(Xax) min(Xax)*1.015 max(Xax)])
ax = gca;
ax.TickDir = 'out';
ax.YDir = 'normal';
xlabel(['\fontsize{13}Track position (cm)' newline '\fontsize{16}Block 2 - Context B'])

suptitle([sDataFiles{fInd}.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sDataFiles{fInd}.imdata(fov).fovLocation newline ...
    'Population vector cross correlation of trial blocks']);

% suptitle('\fontsize{16}Population vector cross correlation of trial blocks');

saveas(gcf,strcat(fullfile(sDataDir,[sDataFiles{fInd}.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sDataFiles{fInd}.imdata(fov).fovLocation '-crossCorrDff']),'.png'));

%}


for fovLoc = 1:1:5
%fovLocation
smoothSpan = 5;
mymap = lines;


if fovLoc == 1
    ind = find([fovLib.isHPC]);
    nameTag = 'HPC';
end
if fovLoc == 2
    ind = find([fovLib.isRSC]);
    nameTag = 'RSC';
end
if fovLoc == 3
    ind = find([fovLib.isPPC]);
    nameTag = 'PPC';
end
if fovLoc == 4
    ind = find([fovLib.isM2]);
    nameTag = 'M2';
end
if fovLoc == 5
    ind = find([fovLib.isV2]);
    nameTag = 'V2';
end

tempA = permute(nanmean(shiftA1A2(:,shift-1:shift+1,ind),2),[1 3 2]);
tempB = permute(nanmean(shiftB1B2(:,shift-1:shift+1,ind),2),[1 3 2]);
tempAB = permute(nanmean(shiftAB(:,shift-1:shift+1,ind),2),[1 3 2]);

figure('Color','white','Position',[0 0 800 330])
suptitle(['\fontsize{13}' 'Mean corr. coeff. ' nameTag ' (n = ' num2str(numel(ind)) ')' newline]);

subplot(1,2,1)
hold on
rectangle('Position',[50,-0.39,140,1.39],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
%rectangle('Position',[49,-0.4,2,1.4],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
%rectangle('Position',[189,-0.4,2,1.4],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
% rectangle('Position',[rewardZones(1,1)-2,0.95,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(1,2)-2,0.95,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(2,1)-2,0.9,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(2,2)-2,0.9,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
plotMeanStd(Xax,tempA,mymap(1,:),smoothSpan);
plotMeanStd(Xax,tempB,mymap(2,:),smoothSpan);
plotMeanStd(Xax,tempAB,mymap(3,:),smoothSpan);


%% plot(Xax,[Cs1(:,shift) Cs4(:,shift) Cs2(:,shift) Cs3(:,shift) Cs5(:,shift) Cs6(:,shift)])

xlabel('\fontsize{13}Track position (cm)')
ylabel('\fontsize{13}Corr. coefficient at diagonal')
%legend('A1 X A2','B1 X A2','A1 X B2','B1 X B2')
%text(110,1.05,'Reward positions','color',[0 1 1])
text(2,-0.3,'Identical')
text(200,-0.3,'Identical')
text(70,-0.3,'Unique VR context')
ax = gca;
ax.TickDir = 'out';
ylim([-0.4 1])


bins = [1:24 96:125];
tempA = permute(nanmean(shiftA1A2(bins,:,ind),1),[2 3 1]);
tempB = permute(nanmean(shiftB1B2(bins,:,ind),1),[2 3 1]);
tempAB = permute(nanmean(shiftAB(bins,:,ind),1),[2 3 1]);


subplot(1,2,2)
hold on

p1 = plotMeanStd((-shift:shift)*2,tempA,mymap(1,:),smoothSpan);
p2 = plotMeanStd((-shift:shift)*2,tempB,mymap(2,:),smoothSpan);
p3 = plotMeanStd((-shift:shift)*2,tempAB,mymap(3,:),smoothSpan);

%plot((-shift:shift)*2,nanmean(Cs6(bins,:)))
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
legend([p1, p2, p3],'A1 X A2','B1 X B2','A X B')
ax = gca;
ax.TickDir = 'out';


saveas(gcf,strcat(fullfile(filePath,[nameTag '-crossCorrCoefDiagonalDff']),'.png'));


%% Begavior plots

figure('Color','white','Position',[0 0 800 330])
suptitle(['\fontsize{13}' 'Behavior signals ' nameTag ' (n = ' num2str(numel(ind)) ')' newline]);

subplot(1,2,1)
hold on
rectangle('Position',[50,0,140,70],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
%rectangle('Position',[49,-0.4,2,1.4],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
%rectangle('Position',[189,-0.4,2,1.4],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
% rectangle('Position',[rewardZones(1,1)-2,0.95,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(1,2)-2,0.95,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(2,1)-2,0.9,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
% rectangle('Position',[rewardZones(2,2)-2,0.9,4,0.05],'FaceColor',[0 1 1],'EdgeColor','none');
plotMeanStd(Xax,velA(:,ind),mymap(1,:),smoothSpan);
plotMeanStd(Xax,velB(:,ind),mymap(2,:),smoothSpan);
%plotMeanStd(Xax,velA2(:,ind),mymap(3,:),smoothSpan);
%plotMeanStd(Xax,velB2(:,ind),mymap(3,:),smoothSpan);

%% plot(Xax,[Cs1(:,shift) Cs4(:,shift) Cs2(:,shift) Cs3(:,shift) Cs5(:,shift) Cs6(:,shift)])

xlabel('\fontsize{13}Track position (cm)')
ylabel('\fontsize{13}Velocity (cm/s)')
%legend('A1 X A2','B1 X A2','A1 X B2','B1 X B2')
%text(110,1.05,'Reward positions','color',[0 1 1])
text(2,3,'Identical')
text(200,3,'Identical')
text(70,3,'Unique VR context')
ax = gca;
ax.TickDir = 'out';
ylim([-0.4 70])



subplot(1,2,2)
hold on
rectangle('Position',[50,0,140,10],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');

p1 = plotMeanStd(Xax,lickA(:,ind),mymap(1,:),smoothSpan);
p2 = plotMeanStd(Xax,lickB(:,ind),mymap(2,:),smoothSpan);


%plot((-shift:shift)*2,nanmean(Cs6(bins,:)))
text(-40,1.05,'Identical part','color',[0 0 0])
%{
subplot(1,2,2)
hold on
plot((-shift:shift)*2,nanmean(Cs1))
plot((-shift:shift)*2,nanmean(Cs2))
plot((-shift:shift)*2,nanmean(Cs3))
plot((-shift:shift)*2,nanmean(Cs4))
%}


xlabel('\fontsize{13}Track position (cm)')
ylabel('\fontsize{13}Lick rate (1/s)')
legend([p1, p2],'Context A','Context B')
ax = gca;
ax.TickDir = 'out';
ylim([-0.04 10])

saveas(gcf,strcat(fullfile(filePath,[nameTag '-behavSignals']),'.png'));

end




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

