% This function concatenates the ROIs for each trial and auto-correlates them on
% trial to trial basis
%
%
%



%{
rois = sData.imdata.activeROIs;
%}


function sData = trialToTrialSimilarity(sData,sDataDir)

fov = 1;
rois = 1:1:size(sData.imdata.binnedRoisDff,3);
trials = 1:1:110;


figure('Color','white','Position',[0 0 800 350])
subplot(1,2,1)
hold on

trackPos = 1:75;
roiMatrix3D = sData.imdata.binnedRoisDff(trials,trackPos,rois);
s  = size(roiMatrix3D);
concatRoiMatrix = fillmissing(reshape(roiMatrix3D,[s(1), s(2)*s(3)]),'constant',0);

imagesc(corr(concatRoiMatrix'))

line([30 30],[1 110],'color','w'); line([1 110],[30 30],'color','w')
line([50 50],[1 110],'color','w'); line([1 110],[50 50],'color','w')
line([70 70],[1 110],'color','w'); line([1 110],[70 70],'color','w')
line([90 90],[1 110],'color','w'); line([1 110],[90 90],'color','w')

title('Opto ON part (0 - 150 cm)')
xlabel('Trials')
ylabel('Trials')
xlim([0 110]); ylim([0 110]); 

subplot(1,2,2)
hold on

trackPos = 76:125;
roiMatrix3D = sData.imdata.binnedRoisDff(trials,trackPos,rois);
s  = size(roiMatrix3D);
concatRoiMatrix = fillmissing(reshape(roiMatrix3D,[s(1), s(2)*s(3)]),'constant',0);

imagesc(corr(concatRoiMatrix'))

line([30 30],[1 110],'color','w'); line([1 110],[30 30],'color','w')
line([50 50],[1 110],'color','w'); line([1 110],[50 50],'color','w')
line([70 70],[1 110],'color','w'); line([1 110],[70 70],'color','w')
line([90 90],[1 110],'color','w'); line([1 110],[90 90],'color','w')

title('Opto OFF part (152 - 250 cm)')
xlabel('Trials')
ylabel('Trials')
xlim([0 110]); ylim([0 110]); 

suptitle(sData.sessionInfo.sessionID)

%saveas(gcf,[fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '_Fov' num2str(fov)  ' ' sData.imdata(fov).fovLocation '_trialToTrialSimilarity_actROIs' nameTag]),'.png']);

saveas(gcf,[fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '_' sData.imdata(fov).fovLocation '_trialToTrialSimilarity']),'.png']);

end


%{
A = [1 1 1; 2 2 2; 3 3 3];

for i = 1:1:3
    B(:,:,i) = A + i*10;
    
    B(:,:,i)
end
reshape(B,[3 9])
%}



