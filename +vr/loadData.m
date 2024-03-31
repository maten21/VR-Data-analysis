function [sDataFiles,filePath] = loadData(method,fileNames,filePath)
%% Load multiple sData files for analysis
% If method = 'light' clears daqdata and ROI signals for faster operation
% and saving RAM memory.
%
%
if ~exist('method')
    method = '';
end

if nargin < 3

if nargin == 2 && iscell(fileNames)

    filePath = '';
else    

% select single or multiple files to analyze 
[fileNames,filePath,~] = uigetfile('*.mat','','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','on' );
end

if ~iscell(fileNames)
    fileNames = {fileNames};
end
    
fileNumber = size(fileNames,2);

sDataFiles = cell(1,fileNumber);

for f = 1:1:fileNumber
    
    fileName = fileNames{f};
    
    load(fullfile(filePath,fileName));
    % insert analysis function here:
    if strcmp(method,'light')
        sData = rmfield(sData,'daqdata');
    if isfield(sData,'imdata') 
    sData.imdata = rmfield(sData.imdata,'roiSignals');
    end
    end
    
    sDataFiles{f} = sData;
    
    clear('sData');
    
    
end




end