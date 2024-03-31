%% Add sessionInfo and mouse info to already made sData struct files
%
%
%
%
%
%
%% Get sData files

clear

% Select and load sData file
[fileNames,sDataDir,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','on' );

if iscell(fileNames) == 0
    fileNames = {fileNames};
end

mouseFolder = 'C:\Users\Mate Neubrandt\Documents\RECORDINGS\MOUSEINFO';


%% READ LAB BOOK DATA FROM TXT FILE 
for i = 1:1:numel(fileNames)
    
load(fullfile(sDataDir,fileNames{i}));

sessionID = strsplit(fileNames{i},'.mat');
sessionID = sessionID{1};
    


% Read lab book info to labBook 
if isfile([sDataDir sessionID '.txt'])
labBook = fileread([sDataDir sessionID '.txt']);
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



% SESSIONINFO
sData.sessionInfo.sessionID = sessionID;
sData.sessionInfo.date = splitText{1};
sData.sessionInfo.sessionNumber = str2double(sessionID(16:17));
sData.sessionInfo.sessionStartTime = start{1};
sData.sessionInfo.sessionStopTime = stop{1};

sData.sessionInfo.recordedData = struct();

if find(strcmp(splitText,'2P')) + 1 == find(strcmp(splitText,'imaging:')) && strcmp(splitText(find(strcmp(splitText,'2P')) + 2),'YES')
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

% MOUSEINFO


if isfile([mouseFolder '\' sessionID(1:5) '.mat'])
load([mouseFolder '\' sessionID(1:5)]);
else
    msgbox(['1) Make sure if "mouseFolder" variable is defined in the script and refers to the path where the mouseinfo files are stored.',char(10),char(10),... 
        '2) Make sure if the mouse info file for this mouse is filled in manually and saved according to the mouse naming standards.'],'Mouseinfo data is not found.');
    return
end

sData.mouseInfo = mouseInfo;
clear('mouseInfo');

save(fullfile(sDataDir,sessionID),'sData');


clear('sData','i','index','labBook','origWeightIndex','sessionID','splitText','start','stop','weightIndex','weightPercIndex');
end
