


function [] = meanOptoEffectOnPCs(sData,sDataDir) 

fov = 1;

%roiSubset = sData.imdata.placeCells(1).placeCells;
roiSubset = sData.imdata.activeROIs;
Xax = sData.behavior.trialMatrices.plotXAxis; 
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;
opto = nanmean(sData.behavior.trialMatrices.optStimMatrix);
optoSignal = round(nanmean(optStimMatrix)./max(nanmean(optStimMatrix)));
optoStart = Xax(find(diff(optoSignal)==1));
optoEnd = Xax(find(diff(optoSignal)==-1)+1);

for i = 1:1:5
Ymax(i) = max(nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{i}(roiSubset,:)));
Ymin(i) = min(nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{i}(roiSubset,:)));
end

Ymax = max(Ymax);
Ymin = min(Ymin);


figure('Color','white','Position',[0 0 800 350]);
suptitle(sData.sessionInfo.sessionID)
subplot(1,2,1)
hold on
Ypos = Ymax*1.05;
rectangle('Position',[optoStart,Ypos,optoEnd-optoStart,Ymax-Ymin*0.05],'FaceColor',[1 0.1 0.1],'EdgeColor','none');

plot(Xax,nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{1}(roiSubset,:)),'LineWidth',2)
plot(Xax,nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{5}(roiSubset,:)),'LineWidth',2)
ylim([Ymin Ymax*1.1])
xlabel('Track position (cm)')
ylabel('Mean DFF of all control place cells')
title('Familiar')

subplot(1,2,2)
hold on
%Ypos = max(nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{3}(placeCells,:)))*1.1;
rectangle('Position',[optoStart,Ypos,optoEnd-optoStart,Ymax-Ymin*0.05],'FaceColor',[1 0.1 0.1],'EdgeColor','none');

plot(Xax,nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{3}(roiSubset,:)),'LineWidth',2)
plot(Xax,nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{4}(roiSubset,:)),'LineWidth',2)
plot(Xax,nanmean(sData.imdata.avBinnedRois.avBinnedRoisDff{2}(roiSubset,:)),'LineWidth',2)
ylim([Ymin Ymax*1.1])
xlabel('Track position (cm)')
legend('control','opto','opto-first exposure')
title('New')

suptitle(sData.sessionInfo.sessionID)

    saveas(gcf,fullfile(sDataDir,['meanOptoEffectOnPCs' '-Fov' num2str(fov) '-' sData.imdata.fovLocation   '.png'] ));
    close(gcf)
    
    
end

