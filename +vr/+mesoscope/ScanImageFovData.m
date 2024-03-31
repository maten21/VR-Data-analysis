%% ScanImage 
% Reads metadata from ScanImage generated *.tif files and returns a
% structure containing FOV specific metadata. Note, this function does not
% run on files lacking "ROI Group Metadata" e. g. recorded by free version 
% of ScanImage.
%
function [fovData, meta, isError] = ScanImageFovData(filePath)
%% Read metadata
isError = false; 

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
    nFrames = nansen.stack.utility.findNumTiffDirectories(tiffInfo);
catch %if findNumTiffDirectories() is not available use slower method
    nFrames = length(imfinfo(filePath));
end

fovData(nFOVs) = struct();
imageLengthWoFlyToLines = 0;
for i = 1:1:nFOVs
    fovData(i).name = ['fov' num2str(i)];
    fovData(i).filePath = filePath;
    fovData(i).imagingSystem = softwareData.imagingSystem;
    fovData(i).fovResXY = roiGroups.imagingRoiGroup.rois(i).scanfields.pixelResolutionXY;
    fovData(i).fovCenterUmXY = roiGroups.imagingRoiGroup.rois(i).scanfields(end).centerXY *objectiveResolution;
    fovData(i).fovSizeUmXY = roiGroups.imagingRoiGroup.rois(i).scanfields(end).sizeXY *objectiveResolution;
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



%% return all metadata
meta = struct();
meta.tiffInfo = tiffInfo;
meta.softwareData = softwareData;
meta.ScanImageRoiGroups = RoiGroups.RoiGroups;
meta.fovData = fovData;

if length(roiGroups.imagingRoiGroup.rois(i).scanfields) > 1
    fovData = struct();
    meta = struct();
    isError = true; 
end


end
