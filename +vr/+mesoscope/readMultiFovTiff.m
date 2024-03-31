% Returns a stitched multiFOV image ScanImage
% frames is an oprional input. frames = [start, stop]
%
function imStack = readMultiFovTiff(filePath,frames)

fovData = vr.mesoscope.ScanImageFovData(filePath);

if nargin < 2

    if fovData(1).nFrames <= 100

        frames = [1, fovData(1).nFrames];

    else
        maxFramesToLoad = 100;
        frames = [1, maxFramesToLoad];
    end
end

if numel(frames) < 2
    frames = [frames, frames];
end

% calculate large FOV resolution
nFOVs = length(fovData);

XcentersPx = zeros([1, nFOVs]);
YcentersPx = zeros([1, nFOVs]);
fovResX = zeros([1, nFOVs]);
fovResY = zeros([1, nFOVs]);
umPerPxX = zeros([1, nFOVs]);
umPerPxY = zeros([1, nFOVs]);

for i = 1:1:nFOVs
    XcentersPx(i) = fovData(i).fovCenterUmXY(1)/fovData(i).umPerPxXY(1);
    YcentersPx(i) = fovData(i).fovCenterUmXY(2)/fovData(i).umPerPxXY(2);
    fovResX(i) = fovData(i).fovResXY(1);
    fovResY(i) = fovData(i).fovResXY(2);
    umPerPxX(i) = fovData(i).umPerPxXY(1);
    umPerPxY(i) = fovData(i).umPerPxXY(2);
end

[xMin, ixMin] = min(XcentersPx);
[xMax, ixMax] = max(XcentersPx);
[yMin, iyMin] = min(YcentersPx);
[yMax, iyMax] = max(YcentersPx);

xRes = round(xMax - xMin + fovResX(ixMin)/2 + fovResX(ixMax)/2);
yRes = round(yMax - yMin + fovResY(iyMin)/2 + fovResY(iyMax)/2);

% blank image stack
imStack = zeros(yRes, xRes, numel(frames(1):frames(end)));

% convert to positive integers for indexing
XcentersPx = round(XcentersPx - min(XcentersPx));
YcentersPx = round(YcentersPx - min(YcentersPx));

for i = 1:1:nFOVs
    
    x = [fovData(i).splitRawImageX1X2Y1Y2(1), fovData(i).splitRawImageX1X2Y1Y2(2)];
    y = [fovData(i).splitRawImageX1X2Y1Y2(3), fovData(i).splitRawImageX1X2Y1Y2(4)];
    im = tiffreadVolume(filePath,'PixelRegion',{y,x,frames});
    imStack(((1:size(im,1)) + YcentersPx(i)), ((1:size(im,2)) + XcentersPx(i)), :) = im;
end

% Note that this version only works if all FOV have the same resolution
% um/PX! ScanImage can record differentresolution in Y axis.

if mean(umPerPxX) ~= mean(umPerPxY)

    [~, stretchDim] = max([mean(umPerPxX), mean(umPerPxY)]);

    if stretchDim == 1
        xStretch = round(mean(umPerPxX)/mean(umPerPxY) * xRes);
        yStretch = yRes;
    else
        xStretch = xRes;
        yStretch = round(mean(umPerPxY)/mean(umPerPxX) * yRes);
    end

    imStack = imresize(imStack,[yStretch xStretch]);
end
