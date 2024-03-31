
clear;

[sDataFilesCtrl, filePath] = vr.loadData('light');
[sDataFilesOpsin, filePath] = vr.loadData('light');

%{
trialType = 2;
for f = 1:1:size(sDataFilesCtrl,trialType)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

for f = 1:1:size(sDataFilesOpsin,trialType)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
%}

%%%%%%% 
figure

myMap = lines(2);
lineColor = [0.7 0.7 0.7];
lineWidth = 0.25;
axLim = [-41 40 -41 40];
subplot(2,2,1)
hold on

trialType = 2;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

rectangle('Position',[-40, 0-lineWidth/2, 80, lineWidth],'FaceColor',lineColor,'EdgeColor','none');
rectangle('Position',[0-lineWidth/2, -40, lineWidth, 80],'FaceColor',lineColor,'EdgeColor','none');


plot(xCtrl,yCtrl,'o')
plot(xOpsin,yOpsin,'o')
plot(mean(xCtrl),mean(yCtrl),'*','color',myMap(1,:))
plot(mean(xOpsin),mean(yOpsin),'*','color',myMap(2,:))

axis(axLim)
xlabel('Velocity shift (cm)')
ylabel('Lick frequency shift (cm)')
title('Full trial')


subplot(2,2,2)
hold on

trialType = 3;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

rectangle('Position',[-40, 0-lineWidth/2, 80, lineWidth],'FaceColor',lineColor,'EdgeColor','none');
rectangle('Position',[0-lineWidth/2, -40, lineWidth, 80],'FaceColor',lineColor,'EdgeColor','none');

plot(xCtrl,yCtrl,'o')
plot(xOpsin,yOpsin,'o')
plot(mean(xCtrl),mean(yCtrl),'*','color',myMap(1,:))
plot(mean(xOpsin),mean(yOpsin),'*','color',myMap(2,:))

axis(axLim)
xlabel('Velocity shift (cm)')
ylabel('Lick frequency shift (cm)')
title('0 cm - 40 cm')


subplot(2,2,3)
hold on

trialType = 4;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

rectangle('Position',[-40, 0-lineWidth/2, 80, lineWidth],'FaceColor',lineColor,'EdgeColor','none');
rectangle('Position',[0-lineWidth/2, -40, lineWidth, 80],'FaceColor',lineColor,'EdgeColor','none');

plot(xCtrl,yCtrl,'o')
plot(xOpsin,yOpsin,'o')
plot(mean(xCtrl),mean(yCtrl),'*','color',myMap(1,:))
plot(mean(xOpsin),mean(yOpsin),'*','color',myMap(2,:))

axis(axLim)
xlabel('Velocity shift (cm)')
ylabel('Lick frequency shift (cm)')
title('80 cm - 120 cm')


subplot(2,2,4)
hold on

trialType = 5;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

rectangle('Position',[-40, 0-lineWidth/2, 80, lineWidth],'FaceColor',lineColor,'EdgeColor','none');
rectangle('Position',[0-lineWidth/2, -40, lineWidth, 80],'FaceColor',lineColor,'EdgeColor','none');

plot(xCtrl,yCtrl,'o')
plot(xOpsin,yOpsin,'o')
plot(mean(xCtrl),mean(yCtrl),'*','color',myMap(1,:))
plot(mean(xOpsin),mean(yOpsin),'*','color',myMap(2,:))

axis(axLim)
xlabel('Velocity shift (cm)')
ylabel('Lick frequency shift (cm)')
title('160 cm - reward')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
axLim = [0 6 -50 20];
shift = 0.05;
xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})

subplot(2,2,1)
hold on

trialType = 2;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

plotBeeswarm(1,xCtrl',myMap(1,:),shift)
plotBeeswarm(2,xOpsin',myMap(2,:),shift)
plot(1,mean(xCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(xOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,yCtrl',myMap(1,:),shift)
plotBeeswarm(5,yOpsin',myMap(2,:),shift)
plot(4,mean(yCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(yOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('full trial')




subplot(2,2,2)
hold on

trialType = 3;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

plotBeeswarm(1,xCtrl',myMap(1,:),shift)
plotBeeswarm(2,xOpsin',myMap(2,:),shift)
plot(1,mean(xCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(xOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,yCtrl',myMap(1,:),shift)
plotBeeswarm(5,yOpsin',myMap(2,:),shift)
plot(4,mean(yCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(yOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('0 - 40 cm')




subplot(2,2,3)
hold on

trialType = 4;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

plotBeeswarm(1,xCtrl',myMap(1,:),shift)
plotBeeswarm(2,xOpsin',myMap(2,:),shift)
plot(1,mean(xCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(xOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,yCtrl',myMap(1,:),shift)
plotBeeswarm(5,yOpsin',myMap(2,:),shift)
plot(4,mean(yCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(yOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('80 - 120 cm')



subplot(2,2,4)
hold on

trialType = 5;
for f = 1:1:size(sDataFilesCtrl,2)
    xCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).velShift;
    yCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end
for f = 1:1:size(sDataFilesOpsin,2)
    xOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).velShift;
    yOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickFreqShift;
end

plotBeeswarm(1,xCtrl',myMap(1,:),shift)
plotBeeswarm(2,xOpsin',myMap(2,:),shift)
plot(1,mean(xCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(xOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,yCtrl',myMap(1,:),shift)
plotBeeswarm(5,yOpsin',myMap(2,:),shift)
plot(4,mean(yCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(yOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('160 cm - reward')

%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('Color','white','Position',[0 0 900 300])
axLim = [0 7 -50 20];
shift = 0.05;
xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
trackLength = 1.8; % track length / 100 (to convert % to cm)

subplot(1,3,1)
hold on

trialType = 2;
for f = 1:1:size(sDataFilesCtrl,2)
    distrShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
    ctrlTrDistrCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).lickQuartiles(2);
    ctrlTrFirstLickCtrl(f) = sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).firstLicksCm(2);
end
for f = 1:1:size(sDataFilesOpsin,2)
    distrShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
    ctrlTrDistrOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).lickQuartiles(2);
    ctrlTrFirstLickOpsin(f) = sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).firstLicksCm(2);
end

plotBeeswarm(1,distrShiftCtrl',myMap(1,:),shift)
plotBeeswarm(2,distrShiftOpsin',myMap(2,:),shift)
plot(1,mean(distrShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(distrShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,firstShiftCtrl',myMap(1,:),shift)
plotBeeswarm(5,firstShiftOpsin',myMap(2,:),shift)
plot(4,mean(firstShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(firstShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('full trial')




subplot(1,3,2)
hold on

trialType = 3;
for f = 1:1:size(sDataFilesCtrl,2)
    distrShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
end
for f = 1:1:size(sDataFilesOpsin,2)
    distrShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
end

plotBeeswarm(1,distrShiftCtrl',myMap(1,:),shift)
plotBeeswarm(2,distrShiftOpsin',myMap(2,:),shift)
plot(1,mean(distrShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(distrShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,firstShiftCtrl',myMap(1,:),shift)
plotBeeswarm(5,firstShiftOpsin',myMap(2,:),shift)
plot(4,mean(firstShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(firstShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('Home-box')




subplot(1,3,3)
hold on

trialType = 4;
for f = 1:1:size(sDataFilesCtrl,2)
    distrShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
end
for f = 1:1:size(sDataFilesOpsin,2)
    distrShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
end

plotBeeswarm(1,distrShiftCtrl',myMap(1,:),shift)
plotBeeswarm(2,distrShiftOpsin',myMap(2,:),shift)
plot(1,mean(distrShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(distrShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,firstShiftCtrl',myMap(1,:),shift)
plotBeeswarm(5,firstShiftOpsin',myMap(2,:),shift)
plot(4,mean(firstShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(firstShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('On-the-path')


%{
subplot(2,2,4)
hold on

trialType = 5;
for f = 1:1:size(sDataFilesCtrl,2)
    distrShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftCtrl(f) = (sDataFilesCtrl{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesCtrl{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
end
for f = 1:1:size(sDataFilesOpsin,2)
    distrShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).lickQuartiles(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).lickQuartiles(2)) * trackLength;
    firstShiftOpsin(f) = (sDataFilesOpsin{1,f}.trials.trialTypesMeta(trialType).firstLicksCm(2) - sDataFilesOpsin{1,f}.trials.trialTypesMeta(1).firstLicksCm(2)) * trackLength;
end

plotBeeswarm(1,distrShiftCtrl',myMap(1,:),shift)
plotBeeswarm(2,distrShiftOpsin',myMap(2,:),shift)
plot(1,mean(distrShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(2,mean(distrShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)
%subplot(1,2,2)
%hold on
plotBeeswarm(4,firstShiftCtrl',myMap(1,:),shift)
plotBeeswarm(5,firstShiftOpsin',myMap(2,:),shift)
plot(4,mean(firstShiftCtrl),'o','color',myMap(1,:),'markerfacecolor',myMap(1,:),'MarkerSize',10)
plot(5,mean(firstShiftOpsin),'o','color',myMap(2,:),'markerfacecolor',myMap(2,:),'MarkerSize',10)
axis(axLim)

xticks([1 2 5 6])
xticklabels({'ctrl.','opsin','ctrl.','opsin'})
title('')

%%%%%%
%}


figure

plot(ctrlTrDistrOpsin,distrShiftOpsin,'o')

%corrcoef(ctrlTrDistrOpsin,distrShiftOpsin)


