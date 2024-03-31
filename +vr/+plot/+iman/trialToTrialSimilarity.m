%{

clear

% Select and load sData file
[sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
load(fullfile(sDataDir,sessionID));
sessionID = strsplit(sessionID,'.mat');
sessionID = sessionID{1};





PCs1 = sData.trials.contextsMeta(1).placeCells;
PCs2 = sData.trials.contextsMeta(2).placeCells;
PCs3 = sData.trials.contextsMeta(3).placeCells;
PCs4 = sData.trials.contextsMeta(4).placeCells;
PCs = union(union(PCs1,PCs3),union(PCs2,PCs4));

nonPCs = setdiff(1:sData.imdata.nROIs,PCs);
PCs2only = setdiff(PCs2,PCs1);

binnedRoisDff = permute(sData.imdata.binnedRoisDff, [1 2 3]);

trackPos = [1:1:15, 90:1:125]; % identical bins to include



%}
function sData = trialToTrialSimilarity(sData,sDataDir)

if sData.behavior.trialMatrices.meta.binNumber > 125
    trackPos = 91:1:160; % identical bins to include
    %trackPos = 1:1:125; % identical bins to include
    nameTag = '180to320cm';    
else
    trackPos = [1:1:15, 90:1:125]; % identical bins to include
    %trackPos = 1:1:125; % identical bins to include
    nameTag = '180to30cm';
end

for fov = 1:1:length(sData.imdata)
rois = sData.imdata(fov).activeROIs;

binnedRoisDff = permute(sData.imdata(fov).binnedRoisDff, [2 1 3]);
binnedRoisDeconv = permute(sData.imdata(fov).binnedRoisDeconv, [2 1 3]);
binnedRoisDeconvRate = permute(sData.imdata(fov).binnedRoisDeconvRate, [2 1 3]);

corrDff = [];
corrDeconv = [];
corrDeconvRate = [];

for r = 1:1:sData.imdata(fov).nROIs
    
    M = fillmissing(binnedRoisDff(:,:,r),'linear'); % fillmissing function works in the vertical dimension by default
    corrDff(:,:,r) = corr(M(trackPos,:));
    
    M = fillmissing(binnedRoisDeconv(:,:,r),'linear'); % fillmissing function works in the vertical dimension by default
    corrDeconv(:,:,r) = corr(M(trackPos,:));

    M = fillmissing(binnedRoisDeconvRate(:,:,r),'linear'); % fillmissing function works in the vertical dimension by default
    corrDeconvRate(:,:,r) = corr(M(trackPos,:));
end

%meanCorr = mean(corrM,3);


Xax = 1:1:size(sData.imdata(fov).binnedRoisDff,1);

figure('Color','white','Position',[0 0 800 800])

subplot(2,2,1)

imagesc(Xax,Xax,nanmean(corrDff(:,:,rois),3))
caxis([-0.1 0.5]);

ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' Dff'])
title('Dff')

subplot(2,2,3)

imagesc(Xax,Xax,nanmean(corrDeconvRate(:,:,rois),3))
caxis([-0.1 0.5]);

ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' Deconv.'])
title('Deconv.')


subplot(2,2,2)
M = sData.behavior.trialMatrices.binVel';
corrVel = corr(M(trackPos,:));

imagesc(Xax,Xax,corrVel)
caxis([-0.5 1]);
ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
title('Running speed')

subplot(2,2,4)
M = sData.behavior.trialMatrices.lickFreqInBin';

corrLick = corr(M(trackPos,:));

imagesc(Xax,Xax,corrLick)
caxis([-0.5 1]);
ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
title('LickRate')



suptitle(['Trial-to-trial correlation ' sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' - ' nameTag])

saveas(gcf,[fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '_Fov' num2str(fov)  ' ' sData.imdata(fov).fovLocation '_trialToTrialSimilarity_actROIs' nameTag]),'.png']);

sData.imdata(fov).trialToTrialCorr.dff182to320cm = corrDff;
sData.imdata(fov).trialToTrialCorr.deconv182to320cm = corrDeconv;
sData.imdata(fov).trialToTrialCorr.deconvRate182to320cm = corrDeconvRate;
sData.imdata(fov).trialToTrialCorr.velocity182to320cm = corrVel;
sData.imdata(fov).trialToTrialCorr.lickRate182to320cm = corrLick;

end



%%

if sData.behavior.trialMatrices.meta.binNumber > 125
    trackPos = 1:1:90; % identical bins to include
    nameTag = '0to180cm';
else
    %trackPos = [1:1:15, 90:1:125]; % identical bins to include
    trackPos = 16:1:89; % identical bins to include
    nameTag = '30to180cm';
end

for fov = 1:1:length(sData.imdata)
rois = sData.imdata(fov).activeROIs;

binnedRoisDff = permute(sData.imdata(fov).binnedRoisDff, [2 1 3]);
binnedRoisDeconv = permute(sData.imdata(fov).binnedRoisDeconv, [2 1 3]);
binnedRoisDeconvRate = permute(sData.imdata(fov).binnedRoisDeconvRate, [2 1 3]);

corrDff = [];
corrDeconv = [];
corrDeconvRate = [];

for r = 1:1:sData.imdata(fov).nROIs
    
    M = fillmissing(binnedRoisDff(:,:,r),'linear'); % fillmissing function works in the vertical dimension by default
    corrDff(:,:,r) = corr(M(trackPos,:));
    
    M = fillmissing(binnedRoisDeconv(:,:,r),'linear'); % fillmissing function works in the vertical dimension by default
    corrDeconv(:,:,r) = corr(M(trackPos,:));

    M = fillmissing(binnedRoisDeconvRate(:,:,r),'linear'); % fillmissing function works in the vertical dimension by default
    corrDeconvRate(:,:,r) = corr(M(trackPos,:));
end

%meanCorr = mean(corrM,3);


Xax = 1:1:size(sData.imdata(fov).binnedRoisDff,1);

figure('Color','white','Position',[0 0 800 800])

subplot(2,2,1)

imagesc(Xax,Xax,nanmean(corrDff(:,:,rois),3))
caxis([-0.1 0.5]);

ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' Dff'])
title('Dff')

subplot(2,2,3)

imagesc(Xax,Xax,nanmean(corrDeconvRate(:,:,rois),3))
caxis([-0.1 0.5]);

ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' Deconv.'])
title('Deconv.')


subplot(2,2,2)
M = sData.behavior.trialMatrices.binVel';
corrVel = corr(M(trackPos,:));

imagesc(Xax,Xax,corrVel)
caxis([-0.5 1]);
ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
title('Running speed')

subplot(2,2,4)
M = sData.behavior.trialMatrices.lickFreqInBin';

corrLick = corr(M(trackPos,:));

imagesc(Xax,Xax,corrLick)
caxis([-0.5 1]);
ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
%title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])
title('LickRate')



suptitle(['Trial-to-trial correlation ' sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' - ' nameTag])

saveas(gcf,[fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '_Fov' num2str(fov)  ' ' sData.imdata(fov).fovLocation '_trialToTrialSimilarity_actROIs' nameTag]),'.png']);

sData.imdata(fov).trialToTrialCorr.dff0to180cm = corrDff;
sData.imdata(fov).trialToTrialCorr.deconv0to180cm = corrDeconv;
sData.imdata(fov).trialToTrialCorr.deconvRate0to180cm = corrDeconvRate;
sData.imdata(fov).trialToTrialCorr.velocity0to180cm = corrVel;
sData.imdata(fov).trialToTrialCorr.lickRate0to180cm = corrLick;

end
end
%% Save each ROI Corr matrix
%{


for fov = 1:1:length(sData.imdata)

binnedRoisDff = permute(sData.imdata(fov).binnedRoisDff, [2 1 3]);
% binnedRoisDff = permute(sData.imdata(fov).binnedRoisDeconv, [2 1 3]);
% binnedRoisDff = permute(sData.imdata(fov).binnedRoisDeconvRate, [2 1 3]);

corrM = [];

for r = 1:1:sData.imdata(fov).nROIs
    
    M = fillmissing(binnedRoisDff(:,:,r),'linear'); % fillmissing function works in the vertical dimension by default
    
    corrM = corr(M(trackPos,:));
    



%meanCorr = mean(corrM,3);


Xax = 1:1:size(sData.imdata(fov).binnedRoisDff,1);

figure('Color','white','Position',[0 0 400 400])

imagesc(Xax,Xax,corrM)
caxis([-0.1 0.5]);

ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])

subFolder = ['trialToTrialCorrIdent' '_Fov' num2str(fov)];
if ~isfolder(fullfile(sDataDir,subFolder)); mkdir(fullfile(sDataDir,subFolder)); end

saveas(gcf,[fullfile(sDataDir,subFolder,[sessionID(1:17) '_Fov' num2str(fov)  ' ' sData.imdata(fov).fovLocation '_trialToTrialSimilarity_Roi' num2str(r)]),'.png']);

close gcf
end
end


%}





%{


M = sData.behavior.trialMatrices.timeInBin';

M = sData.behavior.trialMatrices.lickFreqInBin';


M = sData.behavior.trialMatrices.binVel';

corrM = corr(M(trackPos,:));

figure('Color','white','Position',[0 0 400 400])
imagesc(Xax,Xax,corrM)

ylabel('\fontsize{14}Trials')
xlabel('\fontsize{14}Trials')
title([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation])





figure
plot(mean(meanCorr))

figure
imagesc(corrM)


figure
imagesc(M)

%}