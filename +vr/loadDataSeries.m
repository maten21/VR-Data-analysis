function [sDataFiles,filePath] = loadDataSeries(method)
%% Load multiple sData files for analysis
% If method = 'light' clears daqdata and ROI signals for faster operation
% and saving RAM memory.
%
%
if ~exist('method')
    method = '';
end
% select single or multiple files to analyze 
[fileNames,filePath,~] = uigetfile('*.mat','','MultiSelect','on' );

if ~iscell(fileNames)
    fileNames = {fileNames};
end
    
fileNumber = size(fileNames,2);

%sDataFiles = struct(1,fileNumber);
sDataFiles = struct();

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
    
    sDataFiles(f).sData = sData;
    
    clear('sData');
    
    
end 
S = sDataFiles(1).sData;
if fileNumber >1
for f = 2:1:fileNumber

S = [S, sDataFiles(f).sData];

end
end

sDataFiles = S;
end