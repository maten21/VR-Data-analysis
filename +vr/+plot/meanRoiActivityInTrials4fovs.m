
clear

% Select and load sData file
[sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
load(fullfile(sDataDir,sessionID));
sessionID = strsplit(sessionID,'.mat');
sessionID = sessionID{1};



%figure
%imagesc(A)

%fovLoc = {'PPC', 'M2(RSC)', 'RSC', 'RSC'};
%fovLoc = {'PPC', 'V2', 'RSC', 'M2'};
%fovLoc = {'RSC', 'V2', 'M2(/RSC)',  'PPC(/V2)'};
%fovLoc = {'V2', 'PPC', 'RSC', 'M2'};

%fovLoc = {'V2', 'RSC', 'M2',  'PPC/V2'};

%fovLoc = {'V2(RSC)', 'PPC', 'RSC', 'M2'};


blockStarts = [sData.trials.contextsMeta.blockStart];
nBlocTrials = [sData.trials.contextsMeta.nTrials];
Xlims = [0 sum([sData.trials.contextsMeta.nTrials])+1];


if length(sData.imdata) == 1
    
fig = figure;

fov = 1;
hold on
A = nanmean(sData.imdata(fov).binnedRoisDff,3);
rectangle('Position',[blockStarts(2),min(nanmean(A,2)),nBlocTrials(2),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[blockStarts(4),min(nanmean(A,2)),nBlocTrials(4),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
plot(smoothdata(nanmean(A,2),'gaussian',5),'LineWidth',2)
xlim(Xlims)
ylim([min(nanmean(A,2))-0.01*min(nanmean(A,2)) max(nanmean(A,2))+0.01*min(nanmean(A,2))])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['FOV: ' num2str(fov) ' ' sData.imdata(fov).fovLocation])

    
suptitle(sData.sessionInfo.sessionID(1:17))

saveas(fig,fullfile(sDataDir,['meanActivityAllROIs', '.png']));
    
else


fig = figure;

fov = 1;
subplot(2,2,fov)
hold on
A = nanmean(sData.imdata(fov).binnedRoisDff,3);
rectangle('Position',[blockStarts(2),min(nanmean(A,2)),nBlocTrials(2),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[blockStarts(4),min(nanmean(A,2)),nBlocTrials(4),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
plot(smoothdata(nanmean(A,2),'gaussian',5),'LineWidth',2)
xlim(Xlims)
ylim([min(nanmean(A,2))-0.01*min(nanmean(A,2)) max(nanmean(A,2))+0.01*min(nanmean(A,2))])
%xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['FOV: ' num2str(fov) ' ' sData.imdata(fov).fovLocation])

fov = 2;
subplot(2,2,fov)
hold on
A = nanmean(sData.imdata(fov).binnedRoisDff,3);
rectangle('Position',[blockStarts(2),min(nanmean(A,2)),nBlocTrials(2),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[blockStarts(4),min(nanmean(A,2)),nBlocTrials(4),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
plot(smoothdata(nanmean(A,2),'gaussian',5),'LineWidth',2)
xlim(Xlims)
ylim([min(nanmean(A,2))-0.01*min(nanmean(A,2)) max(nanmean(A,2))+0.01*min(nanmean(A,2))])
%xlabel('Trial number')
%ylabel('Mean DFF of all ROIs')
title(['FOV: ' num2str(fov) ' ' sData.imdata(fov).fovLocation])

fov = 3;
subplot(2,2,fov)
hold on
A = nanmean(sData.imdata(fov).binnedRoisDff,3);
rectangle('Position',[blockStarts(2),min(nanmean(A,2)),nBlocTrials(2),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[blockStarts(4),min(nanmean(A,2)),nBlocTrials(4),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
plot(smoothdata(nanmean(A,2),'gaussian',5),'LineWidth',2)
xlim(Xlims)
ylim([min(nanmean(A,2))-0.01*min(nanmean(A,2)) max(nanmean(A,2))+0.01*min(nanmean(A,2))])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['FOV: ' num2str(fov) ' ' sData.imdata(fov).fovLocation])

fov = 4;
subplot(2,2,fov)
hold on
A = nanmean(sData.imdata(fov).binnedRoisDff,3);
rectangle('Position',[blockStarts(2),min(nanmean(A,2)),nBlocTrials(2),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[blockStarts(4),min(nanmean(A,2)),nBlocTrials(4),max(nanmean(A,2))-min(nanmean(A,2))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
plot(smoothdata(nanmean(A,2),'gaussian',5),'LineWidth',2)
xlim(Xlims)
ylim([min(nanmean(A,2))-0.01*min(nanmean(A,2)) max(nanmean(A,2))+0.01*min(nanmean(A,2))])
xlabel('Trial number')
%ylabel('Mean DFF of all ROIs')
title(['FOV: ' num2str(fov) ' ' sData.imdata(fov).fovLocation])

suptitle(sData.sessionInfo.sessionID)


saveas(fig,fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '-meanActivityAllROIs', '.png']));


end



%% Summarize all sessions
clear
% sDataFiles = vr.loadData;



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

indexes = find([fovLib.isRSC]);
rscData = nan(numel(indexes), 140);
j=1;
for i = indexes

A = nanmean(nanmean(sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).binnedRoisDff,3),2);
s1 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(1).blockStart;
s2 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(2).blockStart;
s3 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(3).blockStart;
s4 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockStart;
s5 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockEnd;

A2 = A(s1:s2-1);
rscData(j,1:numel(A2)) = A2;

A2 = A(s2-10:s3-1);
rscData(j,30:29+numel(A2)) = A2;

A2 = A(s3-10:s4-1);
rscData(j,60:59+numel(A2)) = A2;

A2 = A(s4-10:s5);
rscData(j,90:89+numel(A2)) = A2;

j = j + 1;
end


indexes = find([fovLib.isPPC]);
ppcData = nan(numel(indexes), 140);
j=1;
for i = indexes

A = nanmean(nanmean(sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).binnedRoisDff,3),2);
s1 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(1).blockStart;
s2 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(2).blockStart;
s3 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(3).blockStart;
s4 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockStart;
s5 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockEnd;

A2 = A(s1:s2-1);
ppcData(j,1:numel(A2)) = A2;

A2 = A(s2-10:s3-1);
ppcData(j,30:29+numel(A2)) = A2;

A2 = A(s3-10:s4-1);
ppcData(j,60:59+numel(A2)) = A2;

A2 = A(s4-10:s5);
ppcData(j,90:89+numel(A2)) = A2;

% sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).fovLocation

j = j + 1;
end


indexes = find([fovLib.isM2]);
m2Data = nan(numel(indexes), 140);
j=1;
for i = indexes

A = nanmean(nanmean(sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).binnedRoisDff,3),2);
s1 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(1).blockStart;
s2 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(2).blockStart;
s3 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(3).blockStart;
s4 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockStart;
s5 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockEnd;

A2 = A(s1:s2-1);
m2Data(j,1:numel(A2)) = A2;

A2 = A(s2-10:s3-1);
m2Data(j,30:29+numel(A2)) = A2;

A2 = A(s3-10:s4-1);
m2Data(j,60:59+numel(A2)) = A2;

A2 = A(s4-10:s5);
m2Data(j,90:89+numel(A2)) = A2;

 sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).fovLocation

j = j + 1;
end


indexes = find([fovLib.isV2]);
v2Data = nan(numel(indexes), 140);
j=1;
for i = indexes

A = nanmean(nanmean(sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).binnedRoisDff,3),2);
s1 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(1).blockStart;
s2 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(2).blockStart;
s3 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(3).blockStart;
s4 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockStart;
s5 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockEnd;

A2 = A(s1:s2-1);
v2Data(j,1:numel(A2)) = A2;

A2 = A(s2-10:s3-1);
v2Data(j,30:29+numel(A2)) = A2;

A2 = A(s3-10:s4-1);
v2Data(j,60:59+numel(A2)) = A2;

A2 = A(s4-10:s5);
v2Data(j,90:89+numel(A2)) = A2;

 sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).fovLocation

j = j + 1;
end



indexes = find([fovLib.isHPC]);
hpcData = nan(numel(indexes), 140);
j=1;
for i = indexes

A = nanmean(nanmean(sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).binnedRoisDff,3),2);
s1 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(1).blockStart;
s2 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(2).blockStart;
s3 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(3).blockStart;
s4 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockStart;
s5 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockEnd;

A2 = A(s1:s2-1);
hpcData(j,1:numel(A2)) = A2;

A2 = A(s2-10:s3-1);
hpcData(j,30:29+numel(A2)) = A2;

A2 = A(s3-10:s4-1);
hpcData(j,60:59+numel(A2)) = A2;

A2 = A(s4-10:s5);
hpcData(j,90:89+numel(A2)) = A2;

 sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).fovLocation

j = j + 1;
end


data(1,:) = nanmean(rscData)-min(nanmean(rscData));
data(2,:) = nanmean(ppcData)-min(nanmean(ppcData));
data(3,:) = nanmean(m2Data)-min(nanmean(m2Data));
data(4,:) = nanmean(v2Data)-min(nanmean(v2Data));
data(5,:) = nanmean(hpcData)-min(nanmean(hpcData));


%%



fig = figure;


% subplot(2,3,1)
hold on

rectangle('Position',[40,min(data(:)),30,max(data(:))-min(data(:))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[100,min(data(:)),30,max(data(:))-min(data(:))],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
plot(smoothdata(data','gaussian',5),'LineWidth',2)
xlim([-1 130])
ylim([0-max(data(:))*0.01 max(data(:))])
%xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title([])
legend({'RSC','PPC','M2','V2','HPC'})



fig = figure('Color','white','Position',[0 0 900 600]);

mymap = lines;

subplot(2,3,1)
hold on

fovData = rscData;
range = [min(nanmean(fovData))-(max(nanmean(fovData))-min(nanmean(fovData)))/2, max(nanmean(fovData))+(max(nanmean(fovData))-min(nanmean(fovData)))/2];
rectangle('Position',[40,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[100,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
plotMeanStd(1:1:140,fovData,mymap(1,:),5,2)
xlim([-1 130])
ylim([range(1)-(range(2)-range(1))*0.01 range(2)+(range(2)-range(1))*0.01])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['RSC (n = ' num2str(size(fovData,1)) ' sessions)'])



subplot(2,3,2)
hold on

fovData = ppcData;
range = [min(nanmean(fovData))-(max(nanmean(fovData))-min(nanmean(fovData)))/2, max(nanmean(fovData))+(max(nanmean(fovData))-min(nanmean(fovData)))/2];
rectangle('Position',[40,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[100,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');

plotMeanStd(1:1:140,fovData,mymap(1,:),5,2)
xlim([-1 130])
ylim([range(1)-(range(2)-range(1))*0.01 range(2)+(range(2)-range(1))*0.01])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['PPC (n = ' num2str(size(fovData,1)) ' sessions)'])

subplot(2,3,3)
hold on

fovData = hpcData;
range = [min(nanmean(fovData))-(max(nanmean(fovData))-min(nanmean(fovData)))/2, max(nanmean(fovData))+(max(nanmean(fovData))-min(nanmean(fovData)))/2];
rectangle('Position',[40,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[100,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');

plotMeanStd(1:1:140,fovData,mymap(1,:),5,2)
xlim([-1 130])
ylim([range(1)-(range(2)-range(1))*0.01 range(2)+(range(2)-range(1))*0.01])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['HPC (n = ' num2str(size(fovData,1)) ' sessions)'])


subplot(2,3,4)
hold on

fovData = m2Data;
range = [min(nanmean(fovData))-(max(nanmean(fovData))-min(nanmean(fovData)))/2, max(nanmean(fovData))+(max(nanmean(fovData))-min(nanmean(fovData)))/2];
rectangle('Position',[40,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[100,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');

plotMeanStd(1:1:140,fovData,mymap(1,:),5,2)
xlim([-1 130])
ylim([range(1)-(range(2)-range(1))*0.01 range(2)+(range(2)-range(1))*0.01])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['M2 (n = ' num2str(size(fovData,1)) ' sessions)'])


subplot(2,3,5)
hold on

fovData = v2Data;
range = [min(nanmean(fovData))-(max(nanmean(fovData))-min(nanmean(fovData)))/2, max(nanmean(fovData))+(max(nanmean(fovData))-min(nanmean(fovData)))/2];
rectangle('Position',[40,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[100,range(1),30,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');

plotMeanStd(1:1:140,fovData,mymap(1,:),5,2)
xlim([-1 130])
ylim([range(1)-(range(2)-range(1))*0.01 range(2)+(range(2)-range(1))*0.01])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['V2 (n = ' num2str(size(fovData,1)) ' sessions)'])

suptitle('Mean activity of all ROIs (DFF)')


saveas(fig,fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '-meanActivityAllROIs', '.png']));




%% 

indexes = find([fovLib.isHPC]);
hpcData = nan(numel(indexes), 140);
j=1;
for i = indexes

A = nanmean(nanmean(sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).binnedRoisDff,3),2);
s1 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(1).blockStart;
s2 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(2).blockStart;
s3 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(3).blockStart;
s4 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockStart;
%e = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(4).blockEnd; % last for 4 blocks
s5 = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(5).blockStart;
e = sDataFiles{fovLib(i).fileInd}.trials.contextsMeta(5).blockEnd;

A2 = A(s1:s2-1); hpcData(j,1:numel(A2)) = A2;

A2 = A(s2:s3-1); hpcData(j,s2:s2-1+numel(A2)) = A2;

A2 = A(s3:s4-1); hpcData(j,s3:s3-1+numel(A2)) = A2;

% A2 = A(s4:e); hpcData(j,s4:s4-1+numel(A2)) = A2; % last for 4 blocks

A2 = A(s4:s5-1); hpcData(j,s4:s4-1+numel(A2)) = A2; % last for 5 blocks

A2 = A(s5:e); hpcData(j,s5:s5-1+numel(A2)) = A2; % last for 5 blocks

 sDataFiles{fovLib(i).fileInd}.imdata(fovLib(i).fovInd).fovLocation;

j = j + 1;
end

fig = figure('Color','white','Position',[0 0 400 300]);

mymap = lines;

%subplot(2,3,3)
hold on

fovData = hpcData;
range = [min(nanmean(fovData))-(max(nanmean(fovData))-min(nanmean(fovData)))/2, max(nanmean(fovData))+(max(nanmean(fovData))-min(nanmean(fovData)))/2];
rectangle('Position',[s2,range(1),s3-s2,range(2)-range(1)],'FaceColor',[1 0.95 0.95],'EdgeColor','none'); 
rectangle('Position',[s3,range(1),s4-s3,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % last for 5 blocks
%rectangle('Position',[s4,range(1),e-s4+1,range(2)-range(1)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % last for 4 blocks
rectangle('Position',[s4,range(1),s5-s4,range(2)-range(1)],'FaceColor',[1 0.95 0.95],'EdgeColor','none'); % last for 5 blocks

plotMeanStd(1:1:140,fovData,mymap(1,:),5,2)
xlim([-1 130])
ylim([range(1)-(range(2)-range(1))*0.01 range(2)+(range(2)-range(1))*0.01])
xlabel('Trial number')
ylabel('Mean DFF of all ROIs')
title(['HPC Arch same reward location (n = ' num2str(size(fovData,1)) ' sessions)'])










