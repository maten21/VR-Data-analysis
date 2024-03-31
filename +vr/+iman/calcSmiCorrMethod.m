
clear

% select single or multiple files to analyze 
[fileNames,filePath,~] = uigetfile('*.mat','','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','on' );


if iscell(fileNames) == 0
    fileNames = {fileNames};
end
    
fileNumber = size(fileNames,2);

tic;
for f = 1:1:fileNumber
    
    fileName = fileNames{f};
    

    load(fullfile(filePath,fileName));
    
    %sData = calculateSMICorr(sData);
    
    
    
    
    if length(sData.stats.sessionAvs(1).plotXAxis) == 125
        pos = [1:15 111:125]; % ident
    else
        pos = 91:145; % ident bins
    end
    
for fov = 1:1:length(sData.imdata)
MA = sData.imdata(fov).binnedRoisDff([sData.trials.contextsMeta(1).trials sData.trials.contextsMeta(3).trials],pos,:);
MB = sData.imdata(fov).binnedRoisDff([sData.trials.contextsMeta(2).trials sData.trials.contextsMeta(4).trials],pos,:);

if numel(find(isnan(MA(1,:,:)))) > 5  % exclude first trial if too many NaNs i.e. imaging started late
    MA = MA(2:end,:,:);
end

[corrCoefA, isSignificantA] = vr.getSMICorrAllRois(MA(1:2:size(MA,1),:,:),MA(2:2:size(MA,1),:,:));
[corrCoefB, isSignificantB] = vr.getSMICorrAllRois(MB(1:2:size(MB,1),:,:),MB(2:2:size(MB,1),:,:));
[corrCoefAB, isSignificantAB] = vr.getSMICorrAllRois(MA,MB);

%{
for roi = 1:1:sData.imdata(fov).nROIs
sData.imdata(fov).roiMeta(roi).pos220to30corrCoefA = corrCoefA(roi);
sData.imdata(fov).roiMeta(roi).pos220to30isSignCorrA = isSignificantA(roi);
sData.imdata(fov).roiMeta(roi).pos220to30corrCoefB = corrCoefB(roi);
sData.imdata(fov).roiMeta(roi).pos220to30isSignCorrB = isSignificantB(roi);
sData.imdata(fov).roiMeta(roi).pos220to30corrCoefAB = corrCoefAB(roi);
sData.imdata(fov).roiMeta(roi).pos220to30isSignCorrAB = isSignificantAB(roi);
end
%}

for roi = 1:1:sData.imdata(fov).nROIs
sData.imdata(fov).roiMeta(roi).identPartCorrCoefA = corrCoefA(roi);
sData.imdata(fov).roiMeta(roi).identPartIsSignCorrA = isSignificantA(roi);
sData.imdata(fov).roiMeta(roi).identPartCorrCoefB = corrCoefB(roi);
sData.imdata(fov).roiMeta(roi).identPartIsSignCorrB = isSignificantB(roi);
sData.imdata(fov).roiMeta(roi).identPartCorrCoefAB = corrCoefAB(roi);
sData.imdata(fov).roiMeta(roi).identPartIsSignCorrAB = isSignificantAB(roi);
end


end


if length(sData.stats.sessionAvs(1).plotXAxis) == 125
    pos = 26:90; %unique part
else
    pos = 1:75;
end

for fov = 1:1:length(sData.imdata)
MA = sData.imdata(fov).binnedRoisDff([sData.trials.contextsMeta(1).trials sData.trials.contextsMeta(3).trials],pos,:);
MB = sData.imdata(fov).binnedRoisDff([sData.trials.contextsMeta(2).trials sData.trials.contextsMeta(4).trials],pos,:);

if numel(find(isnan(MA(1,:,:)))) > 5  % exclude first trial if too many NaNs i.e. imaging started late
    MA = MA(2:end,:,:);
end

[corrCoefA, isSignificantA] = vr.getSMICorrAllRois(MA(1:2:size(MA,1),:,:),MA(2:2:size(MA,1),:,:));
[corrCoefB, isSignificantB] = vr.getSMICorrAllRois(MB(1:2:size(MB,1),:,:),MB(2:2:size(MB,1),:,:));
[corrCoefAB, isSignificantAB] = vr.getSMICorrAllRois(MA,MB);

for roi = 1:1:sData.imdata(fov).nROIs
sData.imdata(fov).roiMeta(roi).uniquePartCorrCoefA = corrCoefA(roi);
sData.imdata(fov).roiMeta(roi).uniquePartIsSignCorrA = isSignificantA(roi);
sData.imdata(fov).roiMeta(roi).uniquePartCorrCoefB = corrCoefB(roi);
sData.imdata(fov).roiMeta(roi).uniquePartIsSignCorrB = isSignificantB(roi);
sData.imdata(fov).roiMeta(roi).uniquePartCorrCoefAB = corrCoefAB(roi);
sData.imdata(fov).roiMeta(roi).uniquePartIsSignCorrAB = isSignificantAB(roi);
end

end
    


    save(fullfile(filePath,fileName),'sData');
    clear('sData');
    

   disp([num2str(f) ' / ' num2str(fileNumber) ': ' fileName])
end
toc;
