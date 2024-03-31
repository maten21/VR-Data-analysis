
% Select text file to read 
[fileName,filePath,~] = uigetfile('*.txt','','C:\Users\Mate Neubrandt\Documents\RECORDINGS');

parts = strsplit(fileName,'.'); % Clear".tdms" from Filename.
sessionID = parts{1};

%% READ LAB BOOK DATA FROM TXT FILE 

% Read lab book info to labBook 
if isfile([filePath sessionID '.txt'])
labBook = fileread([filePath sessionID '.txt']);
else
    msgbox('Lab book file containing meta data is required in .txt format with identical file name as the .tdms file.','Lab book file is not found.')
end

splitText = strsplit(labBook);

% Set indexes to extract data from splitText cell array
    weightIndex = 8;
    origWeightIndex = 15;
    weightPercIndex = 5;

index = find(strcmp(splitText, 'Session'));
start = splitText(index(1)+2);
stop = splitText(index(1)+5);

index = find(strcmp(splitText, 'Weight:'));


%% SESSIONINFO
sData.sessionInfo.sessionID = sessionID;
sData.sessionInfo.date = splitText{1};
sData.sessionInfo.sessionNumber = str2double(sessionID(16:17));
sData.sessionInfo.sessionStartTime = start{1};
sData.sessionInfo.sessionStopTime = stop{1};

sData.sessionInfo.recordedData = struct();
if imaging == 1
   sData.sessionInfo.recordedData = {'2P'};    
end

sData.sessionInfo.mouseWeight = str2double(splitText{weightIndex});
sData.sessionInfo.mouseOriginalWeight = str2double(splitText{origWeightIndex});
sData.sessionInfo.mouseWeightPercent = str2double(splitText{weightPercIndex});
sData.sessionInfo.labBook = labBook; 

% Test if the correct values are extracted
if isnan(sData.sessionInfo.mouseWeight) + isnan(sData.sessionInfo.mouseOriginalWeight) + isnan(sData.sessionInfo.mouseWeightPercent) > 0
    msgbox('Modify indexes to extract correct values from splitText cell array.','Incorrect indexes!')
    return
end

%% MOUSEINFO

mouseFolder = 'C:\Users\Mate Neubrandt\Documents\RECORDINGS\MOUSEINFO';

if isfile([mouseFolder '\' sessionID(1:5) '.mat'])
load([mouseFolder '\' sessionID(1:5)]);
else
    msgbox(['1) Make sure if "mouseFolder" variable is defined in the script and refers to the path where the mouseinfo files are stored.',char(10),char(10),... 
        '2) Make sure if the mouse info file for this mouse is filled in manually and saved according to the mouse naming standards.'],'Mouseinfo data is not found.');
    return
end

sData.mouseInfo = mouseInfo;
clear('mouseInfo');



%{
%% SESSIONINFO
sData.sessionInfo.sessionID                     % char      (REQUIRED) A unique ID for the session. It takes the form : mouseName-YYYYMMDD-SS (where SS is the session number).
sData.sessionInfo.date                          % char      (REQUIRED) yyyy.mm.dd
sData.sessionInfo.sessionNumber                 % double    (REQUIRED) The session number part of your sessionID.
sData.sessionInfo.sessionStartTime              % char      (REQUIRED) hh:mm:ss
sData.sessionInfo.sessionStopTime               % char      (REQUIRED) hh:mm:ss
sData.sessionInfo.recordedData                  % cell      (REQUIRED): Either {'2P','LFP','Patch','PupilVideo','SecurityCamera'}. These are only peripheral recordings for your DAQ, i.e. you do not need to add running wheel etc. 

sData.sessionInfo.mouseWeight                   % double    (REQUIRED, for water deprivation experiments)
sData.sessionInfo.mouseOriginalWeight           % double    (REQUIRED, for water deprivation experiments)
sData.sessionInfo.wasSessionAborted             % logic     (Optional) Was the session aborted before you wanted it to?
sData.sessionInfo.recDayNumber                  % double    (Optional) 1,2,3,...,N 
sData.sessionInfo.experimentName                % char      (Optional) name of experiment
sData.sessionInfo.protocol                      % char      (Optional) Tells where a documented protocol of your experiments are located. Fex Google docs etc.
%}
