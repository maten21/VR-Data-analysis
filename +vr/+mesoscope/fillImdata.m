% Include imaging data to sData file - MESOSCOPE
% Open in sequence: 1) sData file, 2) ROI signals, 3) ROIs


clear


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

FOVPos = [nan nan];

%}


processedDataFolder = 'C:\Users\Mate Neubrandt\UIO Physiology Dropbox Dropbox\Lab Data\Mate Neubrandt\Motion corrected mesoscope data\';




% Select and load sData file
[fileName,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
load(fullfile(sDataDir,fileName));

[~, name, ext] = fileparts(fullfile(sDataDir,fileName));

sessionID = name(1:17);
% Select and load extracted ROI signals
% [fileNames,filePath,~] = uigetfile('*.mat','Select ROI Signals','E:\RECORDINGS\OS2\MOTION CORRECTED','MultiSelect','on' );

% Load ROI signals and metadata

d = dir(fullfile(processedDataFolder,sessionID,'FovManagerData'));

for i = 3:length(d)    
    if numel(strfind(d(i).name,'_FOVs')) > 0
        load(fullfile(d(i).folder,d(i).name));        
    end
end

d = dir(fullfile(processedDataFolder,sessionID,'metadata'));

for i = 3:length(d)    
    if numel(strfind(d(i).name,'meta2P')) > 0
        load(fullfile(d(i).folder,d(i).name));        
    end
end


d = dir(fullfile(processedDataFolder,sessionID));

sData.imdata = struct;
for i = 1:1:length(meta2P.fovData)

for f = 3:length(d) 
     if numel(strfind(d(f).name,'FOV')) + numel(strfind(d(f).name,'fov')) > 0 
        if d(f).name(end) == num2str(i)
            d2 = dir(fullfile(d(f).folder,d(f).name,'roisignals/'));
            for g = 3:length(d2)
                load(fullfile(d2(g).folder,d2(g).name)) 
            end
        end
    end
end

    
sData.imdata(i).meta = meta2P.fovData(i);
sData.imdata(i).meta.edgeRelBregma = fovDb.Windows.fovArray(i).edge;
sData.imdata(i).meta.centerRelBregma = mean(fovDb.Windows.fovArray(i).edge);
try
j = str2num(meta2P.fovData(i).imagingSystem(5));
catch
   disp('Imaging Path could not be identified, check indexing.') 
end
sData.imdata(i).meta.tiffInfo = meta2P.tiffInfo(j);
sData.imdata(i).meta.softwareData = meta2P.softwareData(j);
sData.imdata(i).meta.ScanImageRoiGroups = meta2P.ScanImageRoiGroups(j);
% sData.imdata.binned = struct;


[nROIs, nSamples] = size(dff);
sData.imdata(i).roiSignals(2).ch = 'green';
sData.imdata(i).roiSignals(1).ch = 'red';
sData.imdata(i).nROIs = nROIs;
sData.imdata(i).nSamples = nSamples;


sData.imdata(i).roiSignals(2).roif = single(roisMeanFRaw);
sData.imdata(i).roiSignals(2).npilf = single(npilMediF);
sData.imdata(i).roiSignals(2).dff = single(dff);

if exist('ciaDeconv')
    
    deconvRate = NaN(size(ciaDeconv)); % calculate deconvolved activity rate
    for roi = 1:1:size(ciaDeconv,1)
        deconvRate(roi,:) = vr.ifreqDeconv(ciaDeconv(roi,:),meta2P.fovData(i).scanFrameRate,5);
    end
    
    sData.imdata(i).roiSignals(2).deconv = single(ciaDeconv);   % single    Deconvolved signal
    sData.imdata(i).roiSignals(2).denoised = single(ciaDenois);
    sData.imdata(i).roiSignals(2).actRateDeconv = single(deconvRate);   % single    Deconvolved signal
    
end
 
% sData.imdata.roiSignals(2).spikes = single();              % single    Estimate spikes
% sData.imdata.roiSignals(2).roifSubtractedNpil = single();  % single    ROI fluorescence after subtracting neuropil
% sData.imdata.roiSignals(2).dffSubtractedNpil = single();   % single    Delta F/F0 after subtracting neuropil

% load signal extacion options

sData.imdata(i).signalExtractionOptions = options;
% load deconvolution options (if exists)

if exist('ciaOptions')
    sData.imdata(i).signalExtractionOptions = ciaOptions{1,1};
end

% Include ROIs 
sData.imdata(i).roiArray = roi_arr;
%{
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
%}
end

if ~isfolder(fullfile(sDataDir,name))
    mkdir(fullfile(sDataDir,name));
end
save(fullfile(sDataDir,name,fileName),'sData');

