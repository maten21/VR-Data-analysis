function [] = transferAndRenameScanimageFiles


tic
defPath1 = 'E:\Diesel2P\RENAMINGTEST\TempPath_P1';
defPath2 = 'E:\Diesel2P\RENAMINGTEST\TempPath_P2';
%path = "C:\Users\Mate Neubrandt\Dropbox (UIO Physiology Dropbox)\Lab Data\Mate Neubrandt\Diesel2P test recordings\mouseTest_00004.tif";

targerPath = 'E:\Diesel2P\RENAMINGTEST\Target';

p1 = dir(defPath1);
p2 = dir(defPath2);


%{
fovData = ScanImageFovData(path);
tiffInfo = Tiff(path);

datetime("now")

list = {'uguøjl','ihgkt','ogtdkø','lohg','khfsr','srywaye','sdky'};
[indx,tf] = listdlg('ListString',list)
%}


p = vertcat(p1,p2);
j = 1;

for i = 1:1:length(p)

if ~p(i).isdir

expDate = datetime(p(i).date,'InputFormat','dd-MMM-yyyy HH:mm:ss');

expDateStr = datestr(expDate,'yyyy_mm_dd');
expDateStr2 = datestr(expDate,'yyyymmdd');

fovData = ScanImageFovData(fullfile(p(i).folder, p(i).name));


if fovData(1).scanFrameRate  > 5
    tag = '_XYT';
else
    tag = '_fullFov';
end

pathTag = p(i).folder(end-2:end);

name = p(i).name(1:end-10);
ext = p(i).name(end-3:end);
counter = p(i).name(end-8:end-4);
counterShort = p(i).name(end-5:end-4);

try
    if ~strcmp(name(end-7:end),expDateStr2)
        name = [name '-' expDateStr2];
    end
catch
    name = [name '-' expDateStr2];
end

newFileName = [name, '-', counterShort, pathTag, tag, ext];

fileData(j).day = expDateStr;
fileData(j).fileName = fullfile(p(i).folder, p(i).name);
fileData(j).targetName = fullfile(targerPath, expDateStr, name, newFileName);
%fileData(i).
j = j + 1;

end
end




%% Copy files

for f = 1:1:length(fileData)

    sourceFile = fileData(f).fileName;
    targetFile = fileData(f).targetName;


    if ~exist(fileparts(targetFile), 'dir')
        mkdir(fileparts(targetFile));
    end
    
    movefile(sourceFile,targetFile,'f')
    %copyfile(sourceFile,targetFile);
    %delete(sourceFile);

end

disp('Files have been renamed and copied to the target folder.')


toc

end



function fovData = ScanImageFovData(filePath)
%% Read metadata

tiffInfo = Tiff(filePath);

imageResXY = [tiffInfo.getTag(256); tiffInfo.getTag(257)];

softwareData = tiffInfo.getTag(305);
rows = strsplit(softwareData, '\n')';

for i = 1:1:numel(rows)
    try
        eval([rows{i} ';'])
    catch
    end
end

softwareData = SI;

imageTag = tiffInfo.getTag(315);
RoiGroups = jsondecode(imageTag);
roiGroups = RoiGroups.RoiGroups;

%% Generate fovData struct

objectiveResolution = softwareData.objectiveResolution; 
nFOVs = numel(RoiGroups.RoiGroups.imagingRoiGroup.rois);

% This try/catch can be replaced by the function name integrated in NANSEN 
try
    nFrames = findNumTiffDirectories(tiffInfo);
catch %if findNumTiffDirectories() is not available use slower method
    nFrames = length(imfinfo(filePath));
end

fovData(nFOVs) = struct();
imageLengthWoFlyToLines = 0;
for i = 1:1:nFOVs
    fovData(i).filePath = filePath;

    fovData(i).imagingSystem = softwareData.imagingSystem;
    fovData(i).fovResXY = roiGroups.imagingRoiGroup.rois(i).scanfields.pixelResolutionXY;
    fovData(i).fovCenterUmXY = roiGroups.imagingRoiGroup.rois(i).scanfields.centerXY *objectiveResolution;
    fovData(i).fovSizeUmXY = roiGroups.imagingRoiGroup.rois(i).scanfields.sizeXY *objectiveResolution;
    fovData(i).umPerPxXY = fovData(i).fovSizeUmXY ./ fovData(i).fovResXY;
    fovData(i).nFrames = nFrames;
    fovData(i).scanFrameRate = softwareData.hRoiManager.scanFrameRate;
    fovData(i).scanFramePeriodMs = softwareData.hRoiManager.scanFramePeriod;

    imageLengthWoFlyToLines = imageLengthWoFlyToLines + fovData(i).fovResXY(2);
end

nAllFlyToLines = imageResXY(2)-imageLengthWoFlyToLines;

% Distribute fly back/to lines evenly
n = nFOVs+1;
y2 = 0;

for i = 1:1:nFOVs
    flyToLines(i) = round(nAllFlyToLines/(n-i));
    nAllFlyToLines = nAllFlyToLines - flyToLines(i);
    x2 = fovData(i).fovResXY(1);
    x1 = 1;
    
    y2 = y2 + fovData(i).fovResXY(2) + flyToLines(i);
    y1 = y2 - fovData(i).fovResXY(2) + 1;

    fovData(i).splitRawImageX1X2Y1Y2 = [x1, x2, y1, y2];
    fovData(i).splitRawImageCornersXY= [x1, y2; x2, y2; y1, x2; y1, x1];
end

end

