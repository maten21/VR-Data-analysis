
clear;

[fileNames,filePath,~] = uigetfile('*.tif','','','MultiSelect','on' );
files = fullfile(filePath,fileNames);

if ~iscell(files)
    files = {files};
end

fovData = vr.mesoscope.ScanImageFovData(files{1});
for i = 2:1:length(files)
    tempFovData = vr.mesoscope.ScanImageFovData(files{i});
    fovData = [fovData, tempFovData];
end

% Bregma position in the microscopes coordinate system needs to be
% calibrated for each recordings.
bregmaPosUm = [-1000, -3500];

Xdir = 1;
Ydir = -1;



open("C:\Users\Mate Neubrandt\Documents\MATLAB\fovmanager\+brainmap\+paxinos\dorsal_map.fig");
hFig = gcf; 


prevFile= '';

for i = 1:1:length(fovData)
frames = [1, 10];
imAlpha = 0.8;

if ~isequal(prevFile,fovData(i).filePath) %set LUT for every file only
    im = tiffreadVolume(fovData(i).filePath,'PixelRegion',{[1, inf],[1, inf],frames});
    lutMin = quantile(im(:),0.01);
    lutMax = quantile(im(:),0.99);
    prevFile = fovData(i).filePath;
end


x = [fovData(i).splitRawImageX1X2Y1Y2(1), fovData(i).splitRawImageX1X2Y1Y2(2)];
y = [fovData(i).splitRawImageX1X2Y1Y2(3), fovData(i).splitRawImageX1X2Y1Y2(4)];
im = tiffreadVolume(fovData(i).filePath,'PixelRegion',{y,x,frames});


absPos = [0, 0];
absPos(1) = round(fovData(i).fovCenterUmXY(1))*Xdir;
absPos(2) = round(fovData(i).fovCenterUmXY(2))*Ydir;
absPos = round(fovData(i).fovCenterUmXY') - bregmaPosUm;

halfFovSize = fovData(i).fovSizeUmXY'/2;

figPosX = [absPos(1)-halfFovSize(1), absPos(1)+halfFovSize(1)]/1000*Xdir; %invert X axis
figPosY = [absPos(2)-halfFovSize(2), absPos(2)+halfFovSize(2)]/1000*Ydir; %invert Y axis

if frames(2) > frames(1)
    im = mean(im,3);
end

im = flip(im,1); %invert Y axis

hIm = imshow(im,[lutMin, lutMax]);

hIm.XData = figPosX;
hIm.YData = figPosY;
hIm.AlphaData = imAlpha;
% hIm.Axes.YDir = 'normal';
end

set(gca,'YDir','normal')
