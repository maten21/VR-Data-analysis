% USE THIS SCRIPT FOR CONTEXT CHANGE/REMAPPING TASK VERSION (DATA AFTER 2022 MAY)
% Updated 21.10.2021 detect pupil imaging data
% Updated 25.11.2020 do not execute analysis if only 1 trial is recorded
% Complex first analysis and processing of recorded TDMS files of single sessions.
% Input: TDMS file
% Output: sData structure and "session overview" figs are auto saved.
%{
clear;
% OPEN FILE. INPUT: MATLAB Data file, TDMS file recorded by LABVIEW
[fileName,filePath,~] = uigetfile('*.tdms','','C:\Users\Mate Neubrandt\Documents\RECORDINGS' );

%}
% 
% 
function [] = AnalyseRawData_V11_Remap(fileName,filePath)
%SET FOR EACH PLOT!!!!!!

plotFig ='Y'; % Plot figure? (Y or N)
fillInMouseInfo = true; %boolean

viewDistance = 50;
samplingRate = 3000;
% rewardZoneWidth = 10;

parts = strsplit(fileName,'.'); % Clear".tdms" from Filename.
sessionID = parts{1};
clear('parts');

mouseFolder = 'C:\Users\Mate Neubrandt\Dropbox (UIO Physiology Dropbox)\MateData\RECORDINGS\MOUSEINFO';

my_tdms_struct = TDMS_getStruct(fullfile(filePath,fileName));
%% READ DATA FROM TDMS FILE
sData = struct; %make a table with Tdms data converted to matlab data file

sData.daqdata.meta.fs = samplingRate;
sData.daqdata.meta.samplingRate = samplingRate;

sData.daqdata.unityPosition = my_tdms_struct.Untitled.Position_in_Unity.data;
sData.daqdata.wheelDistance = my_tdms_struct.Untitled.Distance__cm_.data;
sData.daqdata.lickSignal = my_tdms_struct.Untitled.Lick_signal.data;
sData.daqdata.waterValve = my_tdms_struct.Untitled.Valve_open.data;
sData.daqdata.lickSignalAnalog = my_tdms_struct.Untitled.Analog_lick_signal.data;

if isfield(my_tdms_struct.Untitled,'Frame_signal') && mean(my_tdms_struct.Untitled.Frame_signal.data) > 0.05 % Consideres imaging data if frame signal is present in 10% of the recording    
    sData.daqdata.frameIndex = my_tdms_struct.Untitled.Frame_signal.data;
    imaging = 1;
else
    imaging = 0;
end

try
    sData.daqdata.frameIndex2 = my_tdms_struct.Untitled.Frame_signal_2.data;
catch
end

if isfield(my_tdms_struct.Untitled,'EyeCam_frame_signal') && mean(my_tdms_struct.Untitled.EyeCam_frame_signal.data) > 0.05 % Consideres imaging data if frame signal is present in 10% of the recording
    
    sData.daqdata.pupilCamFrameIndex = my_tdms_struct.Untitled.EyeCam_frame_signal.data;
    try
    pupilDataFolder = 'C:\Users\Mate Neubrandt\UIO Physiology Dropbox Dropbox\Lab Data\Mate Neubrandt\Raw recordings\PupilData';
    catch
       pupilDataFolder = uigetdir('C:\','Please locate and select "pupilData" folder'); 
    end
    
    pupilExpFolder = fullfile(pupilDataFolder, [sessionID(7:10), '.', sessionID(11:12), '.', sessionID(13:14), '\']); % experimental day folder
    try
        if isfile(fullfile(pupilExpFolder,['-' sessionID '_pupilData.mat']))
            sData = vr.processPupilData(sData,sessionID,pupilExpFolder);
        end
    catch
        msgbox(['ERROR:' sessionID 'Pupil data processing failed!'])
    end
    
    try
        if isfile(fullfile(filePath,['-' sessionID '_pupilData.mat']))
            sData = vr.processPupilData(sData,sessionID,filePath);
        end
    catch
        msgbox(['ERROR:' sessionID 'Pupil data processing failed!'])
    end
    %load(fullfile(filePath,['-' sessionID '_pupilData']));
    % load(fullfile(filePath,['-' sessionID '_frametimes']));
    
    %sData.behavior.pupildata.diameter = Radius*2;
    %sData.behavior.pupildata.center = CenterRotated;
    
end


sData.daqdata.contextIndicators = my_tdms_struct.Untitled.VR_context.data;
% keep opto signals and delet them at the end if not an optogenetic experiment
%sData.daqdata.laserON = my_tdms_struct.Untitled.Laser_ON.data;
try
sData.daqdata.optoSignal = my_tdms_struct.Untitled.Photo_stimulus_signal.data;
catch
end
if isfield(my_tdms_struct.Untitled,'Stimulus_protocol')
sData.daqdata.stimProtIndicators = my_tdms_struct.Untitled.Stimulus_protocol.data;
opto = sum(sData.daqdata.stimProtIndicators) > 0;
else
    opto = false;
end
clear('my_tdms_struct');


%% READ LAB BOOK DATA FROM TXT FILE 

% Read lab book info to labBook 
if isfile([filePath sessionID '.txt'])
labBook = fileread([filePath sessionID '.txt']);
else
    msgbox('Lab book file containing meta data is required in .txt format with identical file name as the .tdms file.','Lab book file is not found.')
end

% Chesk first if only 1 trial is recorded 
splitText = strsplit(labBook,'All trials: ');
splitText = strsplit(splitText{2});
trialCount = str2double(splitText{1});

% 
splitText = strsplit(labBook);

% Set indexes to extract data from splitText cell array
    weightIndex = 10;
    origWeightIndex = 17;
    weightPercIndex = 7;

index = find(strcmp(splitText, 'Session'));
start = splitText(index(1)+2);
stop = splitText(index(1)+5);



% SESSIONINFO
sData.sessionInfo.sessionID = sessionID;
sData.sessionInfo.date = splitText{1};
sData.sessionInfo.sessionNumber = str2double(sessionID(16:17));
sData.sessionInfo.sessionStartTime = start{1};
sData.sessionInfo.sessionStopTime = stop{1};

sData.sessionInfo.recordedData = struct();
%if imaging == 1
%   sData.sessionInfo.recordedData = {'2P'};    
%end
if find(strcmp(splitText,'2P')) + 1 == find(strcmp(splitText,'imaging:')) && strcmp(splitText(find(strcmp(splitText,'2P')) + 2),'YES')
    sData.sessionInfo.recordedData = {'2P'};
end


sData.sessionInfo.mouseWeight = str2double(splitText{weightIndex});
sData.sessionInfo.mouseOriginalWeight = str2double(splitText{origWeightIndex});
sData.sessionInfo.mouseWeightPercent = str2double(splitText{weightPercIndex});
sData.sessionInfo.labBook = labBook; 

% Test if the correct values are extracted
if isnan(sData.sessionInfo.mouseWeight) + isnan(sData.sessionInfo.mouseOriginalWeight) + isnan(sData.sessionInfo.mouseWeightPercent) > 0
    msgbox(['ERROR: ' sessionID newline 'Modify indexes to extract correct values from splitText cell array.','Incorrect indexes!'])
    return
end


% Extract stimulus prtotocols

splitText = strsplit(labBook,'Optical stimulation protocol:');
splitText = strsplit(splitText{2},'VR context:');
splitText = strsplit(splitText{1},'\r\n');

if numel(splitText)<3 splitText = {'\r\n', '0: 10/10 none' ,'\r\n'}; end % fixing LabeView bug of not saving Optical stimulus array

stimProtocolsLabBook(1:numel(splitText)-2) = struct();




for i = 1:1:numel(splitText)-2

    % trialTypeIndicator
    splitText2 = strsplit(splitText{i+1},': ');
    stimProtocolsLabBook(i).trialTypeIndicator = str2double(splitText2{1});
    % proportion
    splitText2 = strsplit(splitText{i+1});
    splitText2 = strsplit(splitText2{2},'/');    
    stimProtocolsLabBook(i).nTrials = str2double(splitText2{1});
    stimProtocolsLabBook(i).proportion = str2double(splitText2{1}) / str2double(splitText2{2});
    % full protocol
    splitText2 = strsplit(splitText{i+1},['/' splitText2{2} ' ']);
    stimProtocolsLabBook(i).protocol = splitText2{2};
    
    if ~isequal(stimProtocolsLabBook(i).protocol,'none')
    % waveform
    splitText2 = strsplit(stimProtocolsLabBook(i).protocol,' from ');
    stimProtocolsLabBook(i).waveform = splitText2{1};
    % from
    splitText2 = strsplit(splitText2{2},' to ');
    if isequal('rew',splitText2{1})
        stimProtocolsLabBook(i).from = 'reward';
    else
        stimProtocolsLabBook(i).from = str2double(splitText2{1});
    end
    % to
    splitText2 = strsplit(splitText2{2},' int. ');
    if isequal('rew',splitText2{1})
        stimProtocolsLabBook(i).to = 'reward';
    else
        stimProtocolsLabBook(i).to = str2double(splitText2{1});
    end
    % intensity
    splitText2 = strsplit(splitText2{2},' %');
    stimProtocolsLabBook(i).intensity = str2num(splitText2{1});
    
    end    
end


% Extract masking prtotocols

splitText = strsplit(labBook,'VR context:');
splitText = strsplit(splitText{2},'Trial types (stimulus protocols):');
splitText = strsplit(splitText{1},'\r\n');

vrContextsLabBook(1:numel(splitText)-2) = struct();

for i = 1:1:numel(splitText)-2
    % trialTypeIndicator
    splitText2 = strsplit(splitText{i+1},': ');
    vrContextsLabBook(i).contextIndicator = str2double(splitText2{1});
    % nTrials
    splitText2 = strsplit(splitText{i+1});
    splitText2 = strsplit(splitText2{2},'/');    
    vrContextsLabBook(i).nTrials = str2double(splitText2{1});
    % name
    splitText2 = strsplit(splitText{i+1},[splitText2{2} ' ']);
    vrContextsLabBook(i).name = splitText2{2};
        
end


% Extract trial type arrays

splitText = strsplit(labBook,'(stimulus protocols):');
splitText = strsplit(splitText{2},'\r\n');
stimList = splitText{2};
contextList = splitText{4};

optStimTypeArrayLabBook(1:numel(stimList)) = NaN;
for i = 1:1:numel(stimList)
    optStimTypeArrayLabBook(i) = str2double(stimList(i));
end

contextTypeArrayLabBook(1:numel(contextList)) = NaN;
for i = 1:1:numel(contextList)
    contextTypeArrayLabBook(i) = str2double(contextList(i));
end



% MOUSEINFO
if fillInMouseInfo
    
mouseFolder = 'C:\Users\Mate Neubrandt\Dropbox (UIO Physiology Dropbox)\MateData\RECORDINGS\MOUSEINFO';

if isfile([mouseFolder '\' sessionID(1:5) '.mat'])
load([mouseFolder '\' sessionID(1:5)]);
else
    msgbox(['1) Make sure if "mouseFolder" variable is defined in the script and refers to the path where the mouseinfo files are stored.',char(10),char(10),... 
        '2) Make sure if the mouse info file for this mouse is filled in manually and saved according to the mouse naming standards.'],'Mouseinfo data is not found.');
    return
end

sData.mouseInfo = mouseInfo;
clear('mouseInfo');

end


%% ANALYZE

if trialCount > 2 % only execute if more than 1 trial
%% Prepare data for analysis
sampleNumber = numel(sData.daqdata.unityPosition);
trackPos = sData.daqdata.unityPosition;


%trackPos = sData.daqdata.unityPosition;
corridorLength = max(trackPos) - min(trackPos);
corridorLength = round(corridorLength/10)*10;
teleportPoint = corridorLength;
% Correct datapoints overrun teleport point
trackPos(trackPos>teleportPoint) = trackPos(trackPos>teleportPoint)-teleportPoint;

% shift the X axis by 50 to not split the identical corridor part
if corridorLength > 251 % shift X axis?
    trackPos = trackPos - 50;
    trackPos(trackPos<0) = trackPos(trackPos<0) + corridorLength;
    
    % once reached the teleport point it should not jump back
    tic
    if length(find(diff(trackPos) > 10)) > 0        
        A = find(diff(trackPos) > 10);
        B = find(diff(trackPos) < -10);
        for i = 1:1:numel(A)        
        C = B - A(i);        
        trackPos(A(i)+1:A(i)+min(C(C>0))) = 0;        
        end        
    end
    toc
    trackPos(end) = 0;

end

% Correct datapoints overrun teleport point
% trackPos(trackPos > corridorLength) = trackPos(trackPos > corridorLength) - corridorLength;

trialStartIndexes = find(diff(trackPos) < -10) +1; % changed from -100 to capture very first trials
nAllTrials = numel(trialStartIndexes)-1;



% rewardZone = ceil(maxRecPos/10)*10 - rewardZoneWidth;



%trackPos(find(diff(positionFromFirstFullTrial) == 0) + 1) = NaN;

% Interpolate corridorPosition values (low Unity frame rate but high sampling rate) 

for i = 1:1:nAllTrials + 1
    if i <= nAllTrials
        j = trialStartIndexes(i);
        k = trialStartIndexes(i+1)-1;
    else
        j = trialStartIndexes(i);
        k = sampleNumber;
    end
        tempPos = trackPos(j:k);
        tempPos(find(diff(tempPos) == 0) + 1) = NaN;
        trackPos(j:k) = fillmissing(tempPos,'linear');
        clear('tempPos');
end

trackPos(trackPos > corridorLength) = corridorLength;
trackPos(trackPos < 0) = 0;

velocity(2:sampleNumber) = diff(smoothdata(sData.daqdata.wheelDistance,'gaussian',samplingRate/10))*samplingRate;
acceleration(2:sampleNumber) = diff(smoothdata(velocity,'gaussian',samplingRate/4))*samplingRate;
%acceleration(2:sampleNumber) = diff(velocity)*samplingRate;

%% LICK DETECTION ANALYSIS
[lickEvents, sData.behavior.lickCount, sData.behavior.licErrorPercentage, lickFig] = vr.sData.analyseLickData(sData.daqdata.lickSignal,samplingRate,15,15); 
sData.behavior.totalLickCount = sum(lickEvents);

try
saveas(gcf,strcat(fullfile(filePath,[sessionID '_lickFig']),'.png'));
close(gcf)
catch
end


lickIfreq = vr.ifreq(lickEvents,samplingRate);



%% Fill data in matrices according to trials 

% produce matrices for plotting
binSize = 2;

binNumber = corridorLength/binSize; % column number
binVel = NaN(nAllTrials,binNumber);
licksInBin = NaN(nAllTrials,binNumber); % matrix to store lick number will be filled cumulatively
lickFreqInBin = NaN(nAllTrials,binNumber);
rewardInBin = NaN(nAllTrials,binNumber); % matrix to store revard events
trialTypeMatrix = NaN(nAllTrials,binNumber);
contextMatrix = NaN(nAllTrials,binNumber);
optStimMatrix = NaN(nAllTrials,binNumber);
accMatrix = NaN(nAllTrials,binNumber);
timeInBin = NaN(nAllTrials,binNumber); % matrix to store delta time values, will be filled cumulatively
%beforeLicksInBin = NaN(nAllTrials,binNumber); % matrix to store predictive licks(drinking excluded)
%timeInBinWheel = NaN(nAllTrials,binNumber); % matrix to store delta time values, will be filled cumulatively


binnedPosition = discretize(trackPos,teleportPoint - corridorLength:binSize:teleportPoint);
% binnedPositionWheel = discretize(corridorPosWheel,-120:binSize:-120+corridorLength);
% Calculate reward events

rewardEvents(1:sampleNumber) = zeros; 
rewardEvents(find(diff(sData.daqdata.waterValve)==1)+1) = 1;

if isfield(sData.daqdata,'stimProtIndicators')
    stimProtIndicators = sData.daqdata.stimProtIndicators;
else
    stimProtIndicators = zeros(size(binnedPosition));
end

if isfield(sData.daqdata,'optoSignal')
    optoSignal = sData.daqdata.optoSignal;
else
    optoSignal = zeros(size(binnedPosition));
end

if isfield(sData.daqdata,'contextIndicators')
    contextIndicators = sData.daqdata.contextIndicators;
else
    contextIndicators = zeros(size(binnedPosition));
end

% Fill all matrices
for r = 1:1:nAllTrials
tempLick = lickEvents(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempLickIfreq = lickIfreq(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempPos = binnedPosition(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempRew = rewardEvents(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempVel = velocity(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempAcc = acceleration(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempContext = contextIndicators(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempOptStim = optoSignal(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempTrialType = stimProtIndicators(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
tempBeforeLick = tempLick;
tempBeforeLick(min(find(tempRew)):numel(tempBeforeLick)) = zeros; 

    for c = 1:1:binNumber
        timeInBin(r,c) = numel(tempPos(tempPos == c));
        licksInBin(r,c) = nansum(tempLick(find(tempPos == c)));
        lickFreqInBin(r,c) = nanmean(tempLickIfreq(find(tempPos == c)));
        rewardInBin(r,c) = nansum(tempRew(find(tempPos == c)));
        binVel(r,c) = nanmean(tempVel(find(tempPos == c)));
        accMatrix(r,c) = nanmean(tempAcc(find(tempPos == c)));
        contextMatrix(r,c) = nanmean(tempContext(find(tempPos == c)));
        optStimMatrix(r,c) = nanmean(tempOptStim(find(tempPos == c)));
        trialTypeMatrix(r,c) = nanmean(tempTrialType(find(tempPos == c)));
        beforeLicksInBin(r,c) = nansum(tempBeforeLick(find(tempPos == c)));
    end

clear('tempPos','tempLick','tempRew','tempVel','tempOptStim','tempTrialType','tempAcc','tempLickIfreq');
end

licksInBin(find(isnan(binVel))) = NaN;
rewardInBin(find(isnan(binVel))) = NaN;
beforeLicksInBin(find(isnan(binVel))) = NaN; % To calculate hit trials
licksPerCm = lickFreqInBin./binVel;

timeInBin = timeInBin/samplingRate;
%binVel = binSize./timeInBin;

%% Analyze trial types, calculate trial numbers and hit rate 

% Identify trial types
%trialTypesCorridor(1:nAllTrials) = NaN;
%{
trialTypesOpto(1:nAllTrials) = NaN;
for i = 1:1:nAllTrials

    j = trialStartIndexes(i);
    k = trialStartIndexes(i+1)-1;
    m = ceil((j+k)/2);
    %trialTypesCorridor(i) = sData.daqdata.trialType(m);
    trialTypesOpto(i) = sData.daqdata.laserON(m);
end
trialTypesCorridor01 = discretize(trialTypesCorridor,2)-1;
%}

% Generate trial type array for contexts
contextTypeArray = mode(contextMatrix,2);
contextTypeIndexes = setdiff(contextTypeArray,-1);



% Generate trial type array for opto
%optStimTypeArray(1:nAllTrials) = NaN;
%%if opto
%    for i = 1:1:nAllTrials
%        optStimTypeArray(i) = sData.daqdata.stimProtIndicators(ceil(mean(trialStartIndexes(i:i+1))));
%    end
%end
optStimTypeArray = mode(trialTypeMatrix,2);
optStimTypeIndexes = setdiff(optStimTypeArray,-1); % get all the used trial type indicators (all are >= 0)


% Calculate hit trials from RewardInBin matrix. 
% A trial is considered correct if rewawd was presented within the reward zone

% rewardZoneWidth = 0;
% RZstart = binNumber - rewardZoneWidth/2 + 1;
% RZend = binNumber;

% hitTrials = nansum(beforeLicksInBin(:,RZstart:RZend),2); % interpreted as lick in RZ
% hitTrials(hitTrials >= 1) = 1; % Correct if there were multiple reward

% rewardedTrials = nansum(rewardInBin(:,RZstart:RZend),2);  % Correction was needed  because of passive lick protocol durin training
% rewardedTrials(rewardedTrials >= 1) = 1; % Correct if there were multiple reward

% predHitTrials = licksPerCm(:,RZstart-1); % Interpreted as minimum 1 lick event/ 10 cm, calculated from lickrate./velocity matrix
% predHitTrials(predHitTrials > 0.1) = 1;
% predHitTrials(predHitTrials <= 0.1) = 0;

% Calculate trial blocks based on both VR context and opt. stim.
A = contextTypeArray*100 +optStimTypeArray;

blockStarts = [1 find(diff(A))'+1];
blockEnds = [find(diff(A))' numel(A)];

isFixedTrialBlocks = numel(blockStarts) == numel(stimProtocolsLabBook);

%blockStarts = [1 find(diff(contextTypeArray))'+1];
%blockEnds = [find(diff(contextTypeArray))' numel(contextTypeArray)];
%contextsMeta(1,numel(blockStarts)) = struct();

try
contextsMeta(1,numel(vrContextsLabBook)) = struct();
for i = 1:1:numel(vrContextsLabBook)
    
    contextsMeta(i).contextIndicator = mean(contextTypeArray(blockStarts(i):blockEnds(i)));
    contextsMeta(i).name = vrContextsLabBook(i).name;
    contextsMeta(i).blockStart = blockStarts(i);
    contextsMeta(i).blockEnd = blockEnds(i);
    contextsMeta(i).trials = blockStarts(i):1:blockEnds(i);
    contextsMeta(i).nTrials = numel(contextsMeta(i).trials);
    
end
catch
end

if isFixedTrialBlocks % I have kept two processing options: if true trial blocks would be fixed with each block allocated with one context and one stimulus protocol. If false, more combination of contexts and stimulations can happen e.g. in case of shuffled stimuluse arrays.
    
    trialTypesMeta(1,numel(stimProtocolsLabBook)) = struct();
for i = 1:1:numel(stimProtocolsLabBook)

    trialTypesMeta(i).trialTypeIndicator = mean(optStimTypeArray(blockStarts(i):blockEnds(i)));
    trialTypesMeta(i).name = stimProtocolsLabBook(i).protocol;
    trialTypesMeta(i).stimProtocol = stimProtocolsLabBook(i).protocol;
    trialTypesMeta(i).blockStart = blockStarts(i);
    trialTypesMeta(i).blockEnd = blockEnds(i);
    trialTypesMeta(i).trials = blockStarts(i):1:blockEnds(i); 
    trialTypesMeta(i).nTrials = numel(trialTypesMeta(i).trials);
 
end

else
for i = 1:1:numel(optStimTypeIndexes)

    trialTypesMeta(i).trialTypeIndicator = optStimTypeIndexes(i);
    
    for j = 1:1:numel(optStimTypeIndexes)
        if stimProtocolsLabBook(j).trialTypeIndicator == optStimTypeIndexes(i)
            trialTypesMeta(i).name = stimProtocolsLabBook(j).protocol;
            trialTypesMeta(i).stimProtocol = stimProtocolsLabBook(j).protocol;
        end
    end
    
    trialTypesMeta(i).stimProtocol = stimProtocolsLabBook(i).protocol;
    trialTypesMeta(i).trials = find(optStimTypeArray == optStimTypeIndexes(i));
    trialTypesMeta(i).nTrials = numel(optStimTypeArray(optStimTypeArray == optStimTypeIndexes(i)));
 
end
end

if isFixedTrialBlocks
    trialBlocksMeta(1,numel(blockStarts)) = struct();
for i = 1:1:numel(blockStarts)
    
    trialBlocksMeta(i).contextIndicator = mean(contextTypeArray(blockStarts(i):blockEnds(i)));
    trialBlocksMeta(i).optStimProtIndicator = mean(optStimTypeArray(blockStarts(i):blockEnds(i)));
    trialBlocksMeta(i).contextName = vrContextsLabBook(i).name;
    trialBlocksMeta(i).optStimProtocol = stimProtocolsLabBook(i).protocol;
    trialBlocksMeta(i).blockStart = blockStarts(i);
    trialBlocksMeta(i).blockEnd = blockEnds(i);
    trialBlocksMeta(i).trials = blockStarts(i):1:blockEnds(i);
    trialBlocksMeta(i).nTrials = numel(trialBlocksMeta(i).trials);
    
end

end


%{
%% ANALYSE LICK DISTRIBUTION
%trialStartIndexes = find(diff(sData.daqdata.unityPosition) < -100);

lickPositions = lickEvents;
lickPositions(lickPositions>0) = sData.daqdata.unityPosition(find(lickPositions>0));

firstLickPositionPerc(1:nAllTrials) = NaN;
firstLickPositionCm(1:nAllTrials) = NaN;

quantiles(nAllTrials,3) = NaN;
lickCountInTrials(nAllTrials) = NaN;
    for i = 1:1:nAllTrials
        temp = lickPositions(trialStartIndexes(i):(trialStartIndexes(i+1)-1));
        quantiles(i,:) = quantile(temp(temp > 0 & temp<rewardZone),[0.25 0.50 0.75])/rewardZone*100;
        lickCountInTrials(i) = numel(temp(temp > 0 & temp<rewardZone));
        
        tempLickPos = temp(temp > 0);
        if numel(tempLickPos)>0                                     
        firstLickPositionCm(i) = tempLickPos(1);                    
        firstLickPositionPerc(i) = tempLickPos(1)/rewardZone*100;   
        end
        
        clear('temp','tempLickPos');
    end
    
sData.stats.lickDistribution.meta.notes = '0.25, 0.50, 0.75 quantiles are calculated for each trials and for the whole session separately. These are normalized lick positions on the track (0 is the corridor start 100 is the reward). I.e.: 25 50 75 values indicate random distribution. ';    
sData.stats.lickDistribution.quantilesTrial = quantiles;
% nanmean(quantiles);
sData.stats.lickDistribution.quantilesSession = quantile(lickPositions(lickPositions > 0 & lickPositions<rewardZone),[0.25 0.50 0.75])/rewardZone*100;
sData.stats.lickDistribution.firstLickPositionCm = firstLickPositionCm;
sData.stats.lickDistribution.firstLickPositionPerc = firstLickPositionPerc;

% Correct NaNs in hit trials when they only licked in the reward zone
nanHitTrials = intersect(find(isnan(quantiles(:,2))),find(hitTrials));
sData.stats.lickDistribution.quantilesTrial(nanHitTrials,:) = 100;

%}

%% CALCULATE SESSION AVERAGES IN EAXH CONTEXT FOR PLOTS 
smoothSpan = 1; % I use a gentle smoothing instead of larger bins
sessionAvs = struct();
sessionAvs.meta.smoothSpan = smoothSpan;
sessionAvs.plotXAxis = (1:binNumber) * binSize;

for i = 1:1:numel(vrContextsLabBook)
    
    sessionAvs(i).avBinVel = smoothdata(nanmean(binVel(contextsMeta(i).trials,:),1),smoothSpan);
    sessionAvs(i).avLicksInBin = smoothdata(nanmean(licksInBin(contextsMeta(i).trials,:),1),smoothSpan);
    sessionAvs(i).avLickFreqInBin = smoothdata(nanmean(lickFreqInBin(contextsMeta(i).trials,:),1),smoothSpan);
    sessionAvs(i).avOptStimMatrix = smoothdata(nanmean(optStimMatrix(contextsMeta(i).trials,:),1),smoothSpan);
    sessionAvs(i).avTimeInBin = smoothdata(nanmean(timeInBin(contextsMeta(i).trials,:),1),smoothSpan);
     
end


% Medians
smoothSpan = 5; % I use a gentle smoothing instead of larger bins
sessionMedians = struct();
sessionMedians.meta.smoothSpan = smoothSpan;
sessionMedians.plotXAxis = (1:binNumber) * binSize;

for i = 1:1:numel(vrContextsLabBook)
    
    sessionMedians(i).medBinVel = smoothdata(nanmedian(binVel(contextsMeta(i).trials,:),1),smoothSpan);
    sessionMedians(i).medLicksInBin = smoothdata(nanmedian(licksInBin(contextsMeta(i).trials,:),1),smoothSpan);
    sessionMedians(i).medLickFreqInBin = smoothdata(nanmedian(lickFreqInBin(contextsMeta(i).trials,:),1),smoothSpan);
    sessionMedians(i).medOptStimMatrix = smoothdata(nanmedian(optStimMatrix(contextsMeta(i).trials,:),1),smoothSpan);
    sessionMedians(i).medTimeInBin = smoothdata(nanmedian(timeInBin(contextsMeta(i).trials,:),1),smoothSpan);

end


sData.stats.sessionAvs = sessionAvs;
sData.stats.sessionMedians = sessionMedians;
%sData.stats.trialNumbers = trialNumbers;
%sData.behavior.trialNumbers = trialNumbers;

%% sData behavior

sData.behavior.signals.corridorPosition = trackPos;
%sData.behavior.signals.corridorPosWheel = corridorPosWheel;
%sData.behavior.signals.lickPositions = lickPositions;
sData.behavior.signals.lickEvents = lickEvents;
sData.behavior.signals.rewardEvents = rewardEvents;
sData.behavior.signals.velocity = velocity;
sData.behavior.signals.lickIfreq = lickIfreq;


%% SAVE MATRICES
matrices = struct();
matrices.meta.binSize = binSize;
matrices.meta.binNumber = binNumber;
matrices.meta.nAllTrials = nAllTrials;
matrices.binVel = binVel;
matrices.accMatrix = accMatrix;
matrices.licksInBin = licksInBin;
%matrices.beforeLicksInBin = beforeLicksInBin;
matrices.timeInBin = timeInBin;
matrices.rewardInBin = rewardInBin;
matrices.contextMatrix = contextMatrix;
matrices.trialTypeMatrix = trialTypeMatrix;
matrices.optStimMatrix = optStimMatrix;
matrices.lickFreqInBin = lickFreqInBin;
matrices.licksPerCm = licksPerCm;
matrices.plotXAxis = (1:binNumber) * binSize; % - 140;
%matrices.beforeLicksInBin = beforeLicksInBin;

sData.behavior.trialMatrices = matrices;

% sData.behavior.rewardZone = rewardZone;
sData.behavior.viewDistance = viewDistance;
%sData.behavior.sessionLengthMin = sum(sum(timeInBin))/60;
sData.behavior.stimProtocols = stimProtocolsLabBook;
%{
hitRates.hitRateNBC = hitRateNBC;
hitRates.hitRateBCN = hitRateBCN;
hitRates.hitRateLightOffBCN = hitRateLightOffBCN;
hitRates.hitRateLightOnBCN = hitRateLightOnBCN;
hitRates.hitRateAfterLightBCN = hitRateAfterLightBCN;
hitRates.hitRateLightOffWithoutAfterLightBCN = hitRateLightOffWithoutAfterLightBCN;
hitRates.hitRateLightOffNBC = hitRateLightOffNBC;
hitRates.hitRateLightOnNBC = hitRateLightOnNBC;
hitRates.hitRateAfterLightNBC = hitRateAfterLightNBC;
hitRates.hitRateLightOffWithoutAfterLightNBC = hitRateLightOffWithoutAfterLightNBC;

sData.stats.hitRates = hitRates;
%}




%% Save trial types
 
if ~isequal(contextTypeArrayLabBook(1:numel(contextTypeArrayLabBook)-1)',contextTypeArray)
msgbox('The trial type array from the labbook is not equal to the array calculated from the DAQ data!');
end

if isFixedTrialBlocks
    sData.trials.trialBlocksMeta = trialBlocksMeta;
end
sData.trials.trialTypeArrayContext = contextTypeArray;
sData.trials.trialTypeArrayStim = optStimTypeArray; % array calculated from the DAQ data!
sData.trials.contextsMeta = contextsMeta;
sData.trials.trialTypesMeta = trialTypesMeta;
sData.trials.vrContextProtocols = vrContextsLabBook;
sData.trials.stimProtocols = stimProtocolsLabBook;

    
%% Save trial types (Updated 2022.11.14)
%{ 
if ~isequal(contextTypeArrayLabBook(1:numel(contextTypeArrayLabBook)-1)',contextTypeArray)
msgbox('The trial type array from the labbook is not equal to the array calculated from the DAQ data!');
end


sData.trials.trialTypeArrayStim = optStimTypeArray; % array calculated from the DAQ data!
sData.trials.trialTypeArrayContext = contextTypeArrayLabBook(1:numel(contextTypeArrayLabBook)-1);
sData.trials.contextsMeta = contextsMeta;
sData.trials.trialTypesMeta = trialTypesMeta;
sData.trials.vrContextProtocols = vrContextsLabBook;
sData.trials.stimProtocols = stimProtocolsLabBook;
%}
%% Calculate some additional data based on trila types

%{
j = find(sData.stats.sessionAvs(1).plotXAxis == binSize); %60/2 + 1; % first bin of the corridor
k = j + rewardZone/binSize - 1; % last bin before RZ 
for t = 1:1:numel(trialTypeIndexes)

    vr.trialType; % short script to calculate trials and nTrials 
    %sData.trials.trialTypesMeta(t).trials = trials;
    sData.trials.trialTypesMeta(t).lickQuartiles = nanmedian(sData.stats.lickDistribution.quantilesTrial(trials,:));
    sData.stats.lickDistribution.trialTypes(t).quartiles = nanmedian(sData.stats.lickDistribution.quantilesTrial(trials,:));
    
    sData.trials.trialTypesMeta(t).firstLicksCm = quantile(sData.stats.lickDistribution.firstLickPositionCm(trials),[0.25 0.50 0.75]);
    sData.trials.trialTypesMeta(t).firstLicksPerc = quantile(sData.stats.lickDistribution.firstLickPositionPerc(trials),[0.25 0.50 0.75]);
    sData.trials.trialTypesMeta(t).accBeforeRZ = quantile(sData.behavior.trialMatrices.accMatrix(trials,k),[0.25 0.50 0.75]);
end
%}


% sData = vr.landmarkModulationIndex(sData);



end % execution of the whole analysis part

%% CLEAR UNNECESSARY FIELDS

if ~opto % Clear opto stimulus data if no stimulated trials are present
    clear('sData.daqdata.laserON','sData.daqdata.optoSignal');
end

save(fullfile(filePath,sessionID),'sData');


%% PLOT FIGURES

if plotFig == 'Y' && trialCount > 2 

nTrialTypes = numel(sData.trials.contextsMeta);
rows = nTrialTypes;
columns = 3;
nSubplots = (nTrialTypes+1)*3;

% for lick plots
myMap = parula(64);
% myMap(1,:) = [0 0 0];

Xax = sData.stats.sessionMedians(1).plotXAxis;

    
Fig1 = figure('Color','white','Position',[0 0 300*columns 150*rows]); %pos of figure [left bottom width height]

    
for t = 1:1:nTrialTypes
    
    %trials = vr.trialType;
    %nTrials = numel(trials);
    %vr.trialType; % 
    trials = sData.trials.contextsMeta(t).trials;
    nTrials = sData.trials.contextsMeta(t).nTrials;
    protType = sData.trials.contextsMeta(t).name;
        
    subplot(rows,columns,(t-1)*3+1);
    imagesc(Xax,1:nTrials,binVel(trials,:)); %(1:number of bins;1:number of trials)
    c = colorbar;
    colormap(gca,hot);
    c.Label.String = 'Speed (cm/s)';
    c.Label.FontSize = 11;
    c.TickDirection = 'out';
    %caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white
    ax = gca;
    ax.TickDir = 'out';
    ylabel('Trials');
    
    subplot(rows,columns,(t-1)*3+2);
%    imagesc(Xax,1:nTrials,licksInBin(trials,:)) %(1:number of bins;1:number of trials)
    imagesc(Xax,1:nTrials,lickFreqInBin(trials,:)) %(1:number of bins;1:number of trials)
    c = colorbar;
    colormap(gca,myMap);
    c.Label.String = 'Lick frequency (Hz)';
    c.Label.FontSize = 11;
    c.TickDirection = 'out';
%    caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white    
    caxis([0 10]); %set limits for color plot, below 1st black, above 2nd white    
    ax = gca;
    ax.TickDir = 'out';
    title(protType);
    
    
    subplot(rows,columns,(t-1)*3+3);
    imagesc(Xax,1:nTrials,rewardInBin(trials,:)) %(1:number of bins;1:number of trials)
    c = colorbar;
    colormap(gca,myMap);
    c.Label.String = 'Reward position';
    c.Label.FontSize = 11;
    c.TickDirection = 'out';
    caxis([0 1]); %set limits for color plot, below 1st black, above 2nd white    
    ax = gca;
    ax.TickDir = 'out';
    
    
end

saveas(gcf,strcat(fullfile(filePath,sessionID),'_1.png'));
close(gcf)


Fig2 = figure('Color','white','Position',[0 0 900 300]); %pos of figure [left bottom width height]

smoothSpan = 5;

subplot(1,2,1);
hold on

% rectangle('Position',[-150,0.5,150,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
% rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:1:nTrialTypes
    trials = sData.trials.contextsMeta(t).trials;
    nTrials = sData.trials.contextsMeta(t).nTrials;
    protType = sData.trials.contextsMeta(t).name;
    plot(Xax,smoothdata(sessionMedians(t).medBinVel,'gaussian',smoothSpan),'Color',map(t,:)); 
end

% axis([-140 (rewardZone + 10) 0 80]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');

% Generate plot legend cell array
for t = 1:1:nTrialTypes
    trials = sData.trials.contextsMeta(t).trials;
    nTrials = sData.trials.contextsMeta(t).nTrials;
    protType = sData.trials.contextsMeta(t).name;
%    if length(protType) > 1
%        strsplit(protType{1})
%    else
%        protType =protType{1};
%    end
    plotLegend{t} = [protType ' (n = ' num2str(nTrials) ')'];
end

subplot(1,2,2);
hold on

% rectangle('Position',[-150,0.05,150,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
% rectangle('Position',[rewardZone,0.05,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:1:nTrialTypes
    trials = sData.trials.contextsMeta(t).trials;
    nTrials = sData.trials.contextsMeta(t).nTrials;
    %plot(Xax,smoothdata(sessionAvs(t).avLicksInBin,'gaussian',smoothSpan),'Color',map(t,:)); 
    plot(Xax,smoothdata(sessionMedians(t).medLickFreqInBin,'gaussian',smoothSpan),'Color',map(t,:));
end

% axis([-140 (rewardZone + 10) 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
legend(plotLegend,'Location','northwest')

saveas(gcf,strcat(fullfile(filePath,[sessionID '_2']),'.png'));
close(gcf)




end






%{




Fig1 = figure('Color','white','Position',[0 0 1200 750]); %pos of figure [left bottom width height]

% plot dimensions
BCNHeight = (0.5*nBCNTrials/nAllTrials);
NBCHeight = (0.5*nNBCTrials/nAllTrials);
subplotSpacing = 0.05;
relativeWidth = (1 - subplotSpacing*4)/3;

Xax = sessionAvs.plotXAxis;
XaxShort = sessionAvs.plotXAxis(1:(corridorLength-80)/binSize);
binBefRZ = (corridorLength-80)/binSize;

%subplot1:
pos = [subplotSpacing (0.475 + NBCHeight) relativeWidth BCNHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
BCNsubplot = subplot('Position',pos);

imagesc(Xax,1:nBCNTrials,binVel(trialType0,:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Speed (virtual cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';
ylabel('Trials');



%Subplot2:
pos = [subplotSpacing 0.4 relativeWidth NBCHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
NBCsubplot = subplot('Position',pos);
imagesc(Xax,1:nNBCTrials,binVel(trialType1,:)) %(1:number of bins;1:number of trials)

c = colorbar;
colormap(gca,hot);
c.Label.String = 'Speed (virtual cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';
ylabel('Trials');


%subplot2.5:

pos = [subplotSpacing 0.125 relativeWidth*0.79 0.2]; %Specify pos as a four-element vector of the form [left bottom width height]. 
VelocitySubplot = subplot('Position',pos);
hold on

rectangle('Position',[-118,0.5,120,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none');
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelBCN,'-',Xax,sessionAvs.avBinVelNBC,'-'); %imagesc(1:BinNumber,1:nNBCTrials,BinVel(TrialType1,:)) %(1:number of bins;1:number of trials)

axis([-120 (binNumber * binSize - 120) 0 60]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');


% for lick plots
myMap = parula(64);
myMap(1,:) = [0 0 0];
%minLick = min(licksInBin)
%maxLick = max(max(licksInBin));
%[0 0 0; parula(5)];
%licksInBin(isnan(licksInBin)) = -1;

%subplot3:
pos = [(subplotSpacing*2+relativeWidth) (0.475 + NBCHeight) relativeWidth BCNHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
LickBCN = subplot('Position',pos);

imagesc(sessionAvs.plotXAxis,1:nBCNTrials,licksInBin(trialType0,:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';


%Subplot4:
pos = [(subplotSpacing*2+relativeWidth) 0.4 relativeWidth NBCHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
LickNBC = subplot('Position',pos);
imagesc(Xax,1:nNBCTrials,licksInBin(trialType1,:)) %(1:number of bins;1:number of trials)

c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';

%subplot4.5:
pos = [(subplotSpacing*2+relativeWidth) 0.125 relativeWidth*0.79 0.2]; %Specify pos as a four-element vector of the form [left bottom width height]. 
LickBeforeSubplot = subplot('Position',pos);
hold on

rectangle('Position',[-118,0.0025,120,0.6],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
rectangle('Position',[rewardZone,0.0025,20,0.6],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none');
rectangle('Position',[(rewardZone-viewDistance),0.0025,0.5,0.6],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(XaxShort,sessionAvs.avLicksInBinBCN(1:binBefRZ),'-',XaxShort,sessionAvs.avLicksInBinNBC(1:binBefRZ),'-'); %imagesc(1:BinNumber,1:nNBCTrials,BinVel(TrialType1,:)) %(1:number of bins;1:number of trials)

axis([-120 (binNumber * binSize - 120) 0 0.6]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

ylabel('Lick before reward (lick/cm)');
legend(['beaconed (n=' num2str(nBCNTrials) ')'],['non-beaconed (n=' num2str(nNBCTrials) ')'],'Location','northwest')




%subplot5:
rewardcolormap = [1 1 1; 0 0 0.5; 0 0.3 0.9; 0 0.5 1; 0 0.7 0.9];
pos = [(subplotSpacing*3+relativeWidth*2) (0.475 + NBCHeight) relativeWidth BCNHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
RewardBCN = subplot('Position',pos);

imagesc(Xax,1:nBCNTrials,rewardInBin(trialType0,:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,rewardcolormap);
c.Label.String = 'Reward count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.5 4.5]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';



%Subplot6:

pos = [(subplotSpacing*3+relativeWidth*2) 0.4 relativeWidth NBCHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
RewardNBC = subplot('Position',pos);
imagesc(Xax,1:nNBCTrials,rewardInBin(trialType1,:)) %(1:number of bins;1:number of trials)

c = colorbar;
colormap(gca,rewardcolormap);
c.Label.String = 'Reward count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.5 4.5]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';
xticklabels = -120:20:(corridorLength-120);
ylabel('Trials');


%subplot6.5:
pos = [(subplotSpacing*3+relativeWidth*2) 0.125 relativeWidth*0.79 0.2]; %Specify pos as a four-element vector of the form [left bottom width height]. 
LickSubplot = subplot('Position',pos);

%{
rectangle('Position',[0.5,0.5,(120/BinSize),60],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
hold on
rectangle('Position',[((RewardZone+120)/BinSize),0.5,(20/BinSize),60],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none');
%}
%fill([0.5,0,(120/BinSize),60],'FaceColor',[1 0 0],'EdgeColor','none');
semilogy(Xax,sessionAvs.avLicksInBinBCN,'-',Xax,sessionAvs.avLicksInBinNBC,'-'); %imagesc(1:BinNumber,1:nNBCTrials,BinVel(TrialType1,:)) %(1:number of bins;1:number of trials)
axis([(-120+binSize) (binNumber * binSize - 120) 0 3]);

xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Lick count (lick/cm)');

saveas(gcf,strcat(fullfile(filePath,sessionID),'.png'));
close(gcf)


%% PLOT OPTO FIGURES

if opto == 1 % Clear opto stimulus data if no stimulated trials are present

% Optogenetics light-ON/light-OFF trial figure
Fig2 = figure('Color','white','Position',[0 0 1200 500]);


VelocitySubplotBCN = subplot(2,3,1);
hold on

rectangle('Position',[-118,0.5,118,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); %R Z
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelBCNLightOff,Xax,sessionAvs.avBinVelBCNLightOn);

axis([-120 (rewardZone+80) 0 60]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');
text((rewardZone-viewDistance),40,'\leftarrow cue visible on horizon');
text(0,45,'\leftarrow corridor starts');
text(-118,50,'"black box"');
text(rewardZone,50,'RZ');
%legend('muscimol','no injection','saline','Location','northwest')

VelocitySubplotNBC = subplot(2,3,4);
hold on

rectangle('Position',[-118,0.5,118,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelNBCLightOff,Xax,sessionAvs.avBinVelNBCLightOn);
axis([-120 (rewardZone+80) 0 60]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');


LickBeforeSubplotBCN = subplot(2,3,2);
hold on

rectangle('Position',[-118,0.0025,118,0.6],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,20,0.6],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.0025,0.5,0.6],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(XaxShort,sessionAvs.avLicksInBinBCNLightOff(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinBCNLightOn(1:binBefRZ));
axis([-120 (rewardZone+80) 0 0.6]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick before reward (lick/cm)');
legend(['light-off (n=' num2str(nLightOffTrialsBCN) ')'],['light-on (n=' num2str(nLightOnTrialsBCN) ')'],'Location','northwest')
title('Beaconed trials')

LickBeforeSubplotNBC = subplot(2,3,5);
hold on

rectangle('Position',[-118,0.0025,118,0.6],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,20,0.6],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.0025,0.5,0.6],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(XaxShort,sessionAvs.avLicksInBinNBCLightOff(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinNBCLightOn(1:binBefRZ));
axis([-120 (rewardZone+80) 0 0.6]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick before reward (lick/cm)');
legend(['light-off (n=' num2str(nLightOffTrialsNBC) ')'],['light-on (n=' num2str(nLightOnTrialsNBC) ')'],'Location','northwest')
title('Non-beaconed trials');


subplot(2,3,3);
imagesc(Xax,1:nAllTrials,optStimMatrix) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Stimulus (V)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white
%title('Light on non-beaconed trials');



subplot(2,3,6);

xticks = 1:1:2;
xticklabels = {'light OFF', 'light ON'};

plot(1:2,[sData.stats.hitRates.hitRateLightOffNBC sData.stats.hitRates.hitRateLightOnNBC],'o-')
axis([0.5 2.5 0 110]);
%set limits for color plot, below 1st black, above 2nd white
%title('Light on non-beaconed trials');
ylabel('Hit rate (%)')
set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);




%{
AllLickSubplotBCN = subplot(2,3,3);
semilogy(Xax,sessionAvs.avLicksInBinBCNLightOff,Xax,sessionAvs.avLicksInBinBCNLightOn);
axis([-120 (rewardZone+80) 0.001 10]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick count (lick/cm)');
legend(['light-off (n=' num2str(nLightOffTrialsBCN) ')'],['light-on (n=' num2str(nLightOnTrialsBCN) ')'],'Location','northwest')


AllLickSubplotNBC = subplot(2,3,6);
semilogy(Xax,sessionAvs.avLicksInBinNBCLightOff,Xax,sessionAvs.avLicksInBinNBCLightOn);
axis([-120 (rewardZone+80) 0.001 10]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick count (lick/cm)');
legend(['light-off (n=' num2str(nLightOffTrialsNBC) ')'],['light-on (n=' num2str(nLightOnTrialsNBC) ')'],'Location','northwest')
%}

saveas(gcf,strcat(fullfile(filePath,[sessionID '_2']),'.png'));
close(gcf)





% Fig 3
Fig3 = figure('Color','white','Position',[0 0 1200 500]);

Xax = (-118:2:(rewardZone+80));
Fig3VelocitySubplotBCN = subplot(2,3,1);
hold on

rectangle('Position',[-117,0.5,118,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelBCNLightOff,Xax,sessionAvs.avBinVelBCNLightOn,Xax,sessionAvs.avBinVelAfterLightBCN,Xax,sessionAvs.avBinVelBCNLightOffWithoutAfterLight);
axis([-120 (rewardZone+80) 0 60]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');
text((rewardZone-viewDistance),40,'\leftarrow cue visible on horizon');
text(0,45,'\leftarrow corridor starts');
text(-118,50,'"black box"');
text(rewardZone,50,'RZ');
%legend('muscimol','no injection','saline','Location','northwest')

Fig3VelocitySubplotNBC = subplot(2,3,4);
hold on
rectangle('Position',[-118,0.5,118,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelNBCLightOff,Xax,sessionAvs.avBinVelNBCLightOn,Xax,sessionAvs.avBinVelAfterLightNBC,Xax,sessionAvs.avBinVelNBCLightOffWithoutAfterLight);
axis([-120 (rewardZone+80) 0 60]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');


Fig3LickBeforeSubplotBCN = subplot(2,3,2);
hold on

rectangle('Position',[-118,0.0025,118,0.6],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,20,0.6],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.0025,0.5,0.6],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(XaxShort,sessionAvs.avLicksInBinBCNLightOff(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinBCNLightOn(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinAfterLightBCN(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinBCNLightOffWithoutAfterLight(1:binBefRZ));
axis([-120 (rewardZone+80) 0 0.6]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick before reward (lick/cm)');
title('Beaconed trials')

Fig3LickBeforeSubplotNBC = subplot(2,3,5);
hold on

rectangle('Position',[-118,0.0025,118,0.6],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,20,0.6],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.0025,0.5,0.6],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(XaxShort,sessionAvs.avLicksInBinNBCLightOff(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinNBCLightOn(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinAfterLightNBC(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinNBCLightOffWithoutAfterLight(1:binBefRZ));
axis([-120 (rewardZone+80) 0 0.6]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick before reward (lick/cm)');
title('Non-beaconed trials');


Fig3AllLickSubplotBCN = subplot(2,3,3);
semilogy(Xax,sessionAvs.avLicksInBinBCNLightOff,Xax,sessionAvs.avLicksInBinBCNLightOn,Xax,sessionAvs.avLicksInBinAfterLightBCN,Xax,sessionAvs.avLicksInBinBCNLightOffWithoutAfterLight);
axis([-120 (rewardZone+80) 0.001 10]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick count (lick/cm)');
legend(strcat('Light-off trials(n=',num2str(nLightOffTrialsBCN),')'),strcat('Light-on trials (n=',num2str(nLightOnTrialsBCN),')'),strcat('Trials after light stimulus (n=',num2str(nAfterLightTrialsBCN),')'),strcat('Light-off,excl. trials after light stim. (n=',num2str(nLightOffWithoutAfterLightTrialsBCN),')'),'Location','northwest')


Fig3AllLickSubplotNBC = subplot(2,3,6);
semilogy(Xax,sessionAvs.avLicksInBinNBCLightOff,Xax,sessionAvs.avLicksInBinNBCLightOn,Xax,sessionAvs.avLicksInBinAfterLightNBC,Xax,sessionAvs.avLicksInBinNBCLightOffWithoutAfterLight);
axis([-120 (rewardZone+80) 0.001 10]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick count (lick/cm)');
legend(strcat('Light-off trials(n=',num2str(nLightOffTrialsNBC),')'),strcat('Light-on trials (n=',num2str(nLightOnTrialsNBC),')'),strcat('Trials after light stimulus (n=',num2str(nAfterLightTrialsNBC),')'),strcat('Light-off,excl. trials after light stim. (n=',num2str(nLightOffWithoutAfterLightTrialsNBC),')'),'Location','northwest')

saveas(gcf,strcat(fullfile(filePath,[sessionID '_3']),'.png'));
close(gcf)







% Various trial type heat maps
Fig4 = figure('Color','white','Position',[0 0 1200 800]);
% Light Off NBC
subplot(4,3,1)
imagesc(Xax,1:numel(intersect(lightOffTrials,trialType0)),licksInBin(intersect(lightOffTrials,trialType0),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white
title('Light off beaconed trials');

% Light Off NBC WO after light trials
subplot(4,3,2)
imagesc(Xax,1:nLightOffWithoutAfterLightTrialsBCN,licksInBin(setdiff(intersect(lightOffTrials,trialType0),afterLightTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white
title('Light off excl. after light trials');

% Light On trials 
subplot(4,3,3)
imagesc(Xax,1:nLightOnTrialsBCN,licksInBin(intersect(trialType0,lightOnTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white
title('Light on beaconed trials');

% Light Off trials 
subplot(4,3,4)
imagesc(Xax,1:numel(intersect(lightOffTrials,trialType0)),binVel(intersect(lightOffTrials,trialType0),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Velocity (cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white


% Light Off BCN
subplot(4,3,5)
imagesc(Xax,1:nLightOffWithoutAfterLightTrialsBCN,binVel(setdiff(intersect(lightOffTrials,trialType0),afterLightTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Velocity (cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white


% Light On trials 
subplot(4,3,6)
imagesc(Xax,1:nLightOnTrialsBCN,binVel(intersect(trialType0,lightOnTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Velocity (cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white


% NBC trials
% Light Off NBC
subplot(4,3,7)
imagesc(Xax,1:numel(intersect(lightOffTrials,trialType1)),licksInBin(intersect(lightOffTrials,trialType1),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white
title('Light off non-beaconed trials');

% Light Off NBC WO after light trials
subplot(4,3,8)
imagesc(Xax,1:nLightOffWithoutAfterLightTrialsNBC,licksInBin(setdiff(intersect(lightOffTrials,trialType1),afterLightTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white
title('Light off excl. after light trials');

% Light On trials 
subplot(4,3,9)
imagesc(Xax,1:nLightOnTrialsNBC,licksInBin(intersect(trialType1,lightOnTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white
title('Light on non-beaconed trials');

% Light Off trials 
subplot(4,3,10)
imagesc(Xax,1:numel(intersect(lightOffTrials,trialType1)),binVel(intersect(lightOffTrials,trialType1),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Velocity (cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white

% Light Off NBC
subplot(4,3,11)
imagesc(Xax,1:nLightOffWithoutAfterLightTrialsNBC,binVel(setdiff(intersect(lightOffTrials,trialType1),afterLightTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Velocity (cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white


% Light On trials 
subplot(4,3,12)
imagesc(Xax,1:nLightOnTrialsNBC,binVel(intersect(trialType1,lightOnTrials),:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Velocity (cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white
%ax = gca; 

saveas(gcf,strcat(fullfile(filePath,[sessionID '_4']),'.png'));
close(gcf)

end

end
%}
end
%{
%PLOT trial types FIGURE
controlFig = figure('Color','white','Position',[10 10 1600 1000]); %pos of figure [left bottom width height]

% create own heatmap:
mymap = [ 1 0 0; 1 1 0]; %2 colors

% plot dimensions
LightOffHeight = (0.5*nLightOffTrials/nAllTrials);
LightOnHeight = (0.5*nLightOnTrials/nAllTrials);
SubplotSpacing = 0.05;
RelativeWidth = (1 - SubplotSpacing*4)/3;

%subplot1:

pos1 = [SubplotSpacing (0.475 + LightOnHeight) RelativeWidth LightOffHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
BCNsubplot = subplot('Position',pos1);

imagesc(1:BinNumber,1:LaserOffNumber,OptStimMatrix(LaserOffTrials,:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Speed (virtual cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
%caxis([auto auto]); %set limits for color plot, below 1st black, above 2nd white

%xlabel(['Position in unity (virtual cm) ']); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
xticklabels = -120:20:(CorridorLength-120);
%xticklabels = 0:10:(200/binSize); %cannot set the axis to show 157 as the end%
xticks = linspace(1, size(OptStimMatrix, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

ylabel('Trials');
yticklabels = 0:5:LaserOffNumber;
yticks = linspace(1, LaserOffNumber, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)



%Subplot2:

pos1 = [SubplotSpacing 0.4 RelativeWidth LightOnHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
NBCsubplot = subplot('Position',pos1);
imagesc(1:BinNumber,1:TLaserOnNumber,OptStimMatrix(LaserOnTrials,:)) %(1:number of bins;1:number of trials)

c = colorbar;
colormap(gca,hot);
c.Label.String = 'Speed (virtual cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 

%caxis([auto auto]); %set limits for color plot, below 1st black, above 2nd white

%xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
xticklabels = -120:20:(CorridorLength-120);
%xticklabels = 0:10:(200/binSize); %cannot set the axis to show 157 as the end%
xticks = linspace(1, size(OptStimMatrix, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

ylabel('Trials');
yticklabels = 0:5:LaserOnNumber;
yticks = linspace(1, LaserOnNumber, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

%}
