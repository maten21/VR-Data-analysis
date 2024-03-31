function sData = readLabViewMetaData(sData, filePath, sessionID)



% Read lab book info to labBook
if isfile([filePath sessionID '.txt'])
    labBook = fileread([sData '.txt']);
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

stimProtocols(1:numel(splitText)-2) = struct();

for i = 1:1:numel(splitText)-2

    % trialTypeIndicator
    splitText2 = strsplit(splitText{i+1},': ');
    stimProtocols(i).trialTypeIndicator = str2double(splitText2{1});
    % proportion
    splitText2 = strsplit(splitText{i+1});
    splitText2 = strsplit(splitText2{2},'/');    
    stimProtocols(i).proportion = str2double(splitText2{1}) / str2double(splitText2{2});
    % full protocol
    splitText2 = strsplit(splitText{i+1},['/' splitText2{2} ' ']);
    stimProtocols(i).protocol = splitText2{2};
    
    if ~isequal(stimProtocols(i).protocol,'none')
    % waveform
    splitText2 = strsplit(stimProtocols(i).protocol,' from ');
    stimProtocols(i).waveform = splitText2{1};
    % from
    splitText2 = strsplit(splitText2{2},' to ');
    if isequal('rew',splitText2{1})
        stimProtocols(i).from = 'reward';
    else
        stimProtocols(i).from = str2double(splitText2{1});
    end
    % to
    splitText2 = strsplit(splitText2{2},' int. ');
    if isequal('rew',splitText2{1})
        stimProtocols(i).to = 'reward';
    else
        stimProtocols(i).to = str2double(splitText2{1});
    end
    % intensity
    splitText2 = strsplit(splitText2{2},' %');
    stimProtocols(i).intensity = str2num(splitText2{1});
    
    end    
end


% Extract masking prtotocols

splitText = strsplit(labBook,'VR context:');
splitText = strsplit(splitText{2},'Trial types (stimulus protocols):');
splitText = strsplit(splitText{1},'\r\n');

vrContexts(1:numel(splitText)-2) = struct();

for i = 1:1:numel(splitText)-2
    % trialTypeIndicator
    splitText2 = strsplit(splitText{i+1},': ');
    vrContexts(i).contextIndicator = str2double(splitText2{1});
    % proportion
    splitText2 = strsplit(splitText{i+1});
    splitText2 = strsplit(splitText2{2},'/');    
    vrContexts(i).nTrials = str2double(splitText2{1});
    % name
    splitText2 = strsplit(splitText{i+1},[splitText2{2} ' ']);
    vrContexts(i).name = splitText2{2};
        
end


% Extract trial type arrays

splitText = strsplit(labBook,'(stimulus protocols):');
splitText = strsplit(splitText{2},'\r\n');
stimList = splitText{2};
contextList = splitText{4};

trialTypeArrayStim(1:numel(stimList)) = NaN;
for i = 1:1:numel(stimList)
    trialTypeArrayStim(i) = str2double(stimList(i));
end

trialTypeArrayContext(1:numel(contextList)) = NaN;
for i = 1:1:numel(contextList)
    trialTypeArrayContext(i) = str2double(contextList(i));
end
