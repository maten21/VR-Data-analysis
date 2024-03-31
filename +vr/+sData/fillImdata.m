% Include imaging data secorded with OS2 to sData file 
% Open in sequence: 1) sData file, 2) ROI signals, 3) ROIs


% clear

function [] = fillImdata(sDataPath,roiSignalPath,nansenFormat)

FOVPos = [nan nan];
%{ 
FOVPos = [0 0];

% CONVERT OS2 COORDINATES

XCoordOS2 = 3133;
YCoordOS2 = -1407;

windowCenterPos = [0, -2.2]; innerWindowSize = 2.5; % RSC standard
windowCenterPos = [0, 0.2]; innerWindowSize = 2.5; % M2 
windowCenterPos = [0, -1.5]; innerWindowSize = 5; % Large window start at +1

windowAnteriorPos = [windowCenterPos(1), (windowCenterPos(2) + innerWindowSize/2)];

FOVPos = [windowAnteriorPos(1)-YCoordOS2/1000, windowAnteriorPos(2)-XCoordOS2/1000];



%}
if nargin < 3
    nansenFormat = false;
end
% Select and load sData file
if nargin < 2
    sDataPath = '';
    %[sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
    %sDataPath = fullfile(sDataDir,sessionID);
    %load(sDataPath);
    % sessionID = strsplit(sessionID,'.mat');
    % sessionID = sessionID{1};
end
try
    load(sDataPath);    
catch
    [sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
    sDataPath = fullfile(sDataDir,sessionID);
    load(sDataPath);
end

[sDataDir, sessionID, ~] = fileparts(sDataPath);



% Select and load extracted ROI signals
if nargin < 2
    [fileNames,filePath,~] = uigetfile('*.mat',['Select ROI Signals for ' sessionID],'E:\RECORDINGS\OS2\MOTION CORRECTED','MultiSelect','on' );
else
    filePath = fileparts(roiSignalPath);
    d = dir(filePath);
    fileNames = {d(~[d.isdir]).name};    
end
    


if iscell(fileNames) == 0
    fileNames = {fileNames};
end

fileNumber = size(fileNames,2);

for f = 1:1:fileNumber
    [~, ~, ext] = fileparts(fileNames{f});
    if strcmp(ext,'.mat')
        load(fullfile(filePath,fileNames{f}));
    end
end

% move from ROI signals subfolder to the exp folder
filePath = fileparts(fileparts(roiSignalPath));

% Load metadata
files = dir(fullfile(filePath,'metadata'));
load(fullfile(files(3).folder,files(3).name)); 


sData.imdata = struct;
sData.imdata.meta = struct;
% sData.imdata.binned = struct;

try
    s = size(dff);
catch
    s = size(RoiSignals_Dff);
end
nROIs = min(s);
nSamples = max(s);

sData.imdata.roiSignals(2).ch = 'green';
sData.imdata.roiSignals(1).ch = 'red';
sData.imdata.nROIs = nROIs;
sData.imdata.nSamples = nSamples;

try
    if exist('roisMeanFRaw')
        sData.imdata.roiSignals(2).roif = single(roisMeanFRaw);
    elseif exist('roiMeanF')
        sData.imdata.roiSignals(2).roif = single(roiMeanF);
    elseif exist('RoiSignals_MeanF')
        sData.imdata.roiSignals(2).roif = single(RoiSignals_MeanF);
    else
        msgbox('ROI mean signal has not been found.')
    end
catch
    msgbox('ROI mean signal has not been found.')
end
try
    if exist('npilMediF')
    sData.imdata.roiSignals(2).npilf = single(npilMediF);
    elseif exist('neurpilMeanF')
        sData.imdata.roiSignals(2).npilf = single(neurpilMeanF);
    elseif exist('RoiSignals_NeuropilF')
        sData.imdata.roiSignals(2).npilf = single(RoiSignals_NeuropilF);
    else
        msgbox('Neuropil signal has not been found.')
    end
catch
    msgbox('Neuropil signal has not been found.')
end

if exist('dff')
    sData.imdata.roiSignals(2).dff = single(dff);
elseif exist('RoiSignals_Dff')
    sData.imdata.roiSignals(2).dff = single(RoiSignals_Dff);
end
    
if exist('ciaDeconv')
    
    deconvRate = NaN(size(ciaDeconv)); % calculate deconvolved activity rate
    for roi = 1:1:size(ciaDeconv,1)
        deconvRate(roi,:) = vr.ifreqDeconv(ciaDeconv(roi,:),meta2P.fps,10);
    end
    
    sData.imdata.roiSignals(2).deconv = single(ciaDeconv);   % single    Deconvolved signal
    sData.imdata.roiSignals(2).denoised = single(ciaDenois);
    sData.imdata.roiSignals(2).actRateDeconv = single(deconvRate);   % single    Deconvolved signal
    
end

if exist('deconvolved')
    
    deconvRate = NaN(size(deconvolved)); % calculate deconvolved activity rate
    for roi = 1:1:size(deconvolved,1)
        deconvRate(roi,:) = vr.ifreqDeconv(deconvolved(roi,:),meta2P.fps,10);
    end
    
    sData.imdata.roiSignals(2).deconv = single(deconvolved);   % single    Deconvolved signal
    sData.imdata.roiSignals(2).denoised = single(denoised);
    sData.imdata.roiSignals(2).actRateDeconv = single(deconvRate);   % single    Deconvolved signal
    
end

if exist('RoiSignals_Deconvolved')
    
    deconvRate = NaN(size(RoiSignals_Deconvolved)); % calculate deconvolved activity rate
    for roi = 1:1:size(RoiSignals_Deconvolved,1)
        deconvRate(roi,:) = vr.ifreqDeconv(RoiSignals_Deconvolved(roi,:),meta2P.fps,10);
    end
    
    sData.imdata.roiSignals(2).deconv = single(RoiSignals_Deconvolved);   % single    Deconvolved signal
    sData.imdata.roiSignals(2).denoised = single(RoiSignals_Denoised);
    sData.imdata.roiSignals(2).actRateDeconv = single(deconvRate);   % single    Deconvolved signal
    sData.imdata.roiSignals(2).deconvOptions = OptionsDeconvolution;
    
    
end


% sData.imdata.roiSignals(2).spikes = single();              % single    Estimate spikes
% sData.imdata.roiSignals(2).roifSubtractedNpil = single();  % single    ROI fluorescence after subtracting neuropil
% sData.imdata.roiSignals(2).dffSubtractedNpil = single();   % single    Delta F/F0 after subtracting neuropil

% load signal extacion options
if exist('ciaOptions')
sData.imdata.signalExtractionOptions = options;
end
if exist('OptionsSignalExtraction')
sData.imdata.signalExtractionOptions = OptionsSignalExtraction;
end


% load deconvolution options (if exists)

if exist('ciaOptions')
    sData.imdata.signalExtractionOptions = ciaOptions{1,1};
end

try
    % Include ROIs
    sData.imdata.roiArray = roi_arr;
catch
    % Select and load extracted ROI signals
    [fileName,filePath,~] = uigetfile('*.mat',['Select ROI Array for ' sessionID] ,filePath,'MultiSelect','off' );
    load(fullfile(filePath,fileName));
    try
        sData.imdata.roiArray = roi_arr;
    catch
        try
            sData.imdata.roiArray = roiArray;
        catch
            sData.imdata.roiArray = RoiArray;
        end
    end
end

% Manual input of some data
prompt = {'Laser Power (mW)','Wave Length (nm)','FOV Coord X','FOV Coord Y','PMT Gain Green'};
dlgtitle = 'Specify These Imaging Parameters';
definput = {'','930',num2str(FOVPos(1)),num2str(FOVPos(2)),'25'};
answer = inputdlg(prompt,dlgtitle,[1 40],definput);

laserPower = str2double(answer{1});
waveLength = str2double(answer{2});
fovCoordX = str2double(answer{3});
fovCoordY = str2double(answer{4});
fovCoord = [fovCoordX fovCoordY meta2P.zPosition/1000];
pmtGainGreen = str2double(answer{5});


% add metadata
sData.imdata.meta = meta2P;
sData.imdata.meta.laserPower =  laserPower;          % mW                
sData.imdata.meta.waveLength =  waveLength;         % nm          
sData.imdata.meta.fovCoordinates = fovCoord;        % from center of window, AP ML (0 -500 means imaging at -2.2 mm left hemisphere)
sData.imdata.meta.pmtGain(1,1) = pmtGainGreen;


sData = vr.transposeRoiSignals(sData); %transpose signal if it is in different order

if nansenFormat
    filePath = fullfile(sDataDir,['session-' sessionID(1:17)]);
else
    filePath = fullfile(sDataDir,sessionID);
end

if ~exist(filePath, 'dir'); mkdir(filePath);  end

save(fullfile(filePath,sessionID),'sData');

end