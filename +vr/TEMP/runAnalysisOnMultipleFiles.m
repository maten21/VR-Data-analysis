clear
clc

% select single or multiple files to analyze 
[fileNames,filePath,~] = uigetfile('*.tdms','','C:\','MultiSelect','on' );


if iscell(fileNames) == 0
    fileNames = {fileNames};
end
    
fileNumber = size(fileNames,2);

tic;
for f = 1:1:fileNumber
    
    fileName = fileNames{f};
    
    % insert analysis function here:
    if numel(strsplit(fileName,'cm')) > 1
        AnalyseRawData_V10(fileName,filePath);
    else
        vr.sData.AnalyseRawData_V11_Remap(fileName,filePath);
    end

   disp([num2str(f) ' / ' num2str(fileNumber) ': ' fileName])
end
toc;


%% Copy files to target folders

fileStruct = dir(filePath);

targetPaths(1:length(fileStruct)) = {''};
for f = 1:1:length(fileStruct)
    if numel(fileStruct(f).name) > 2 % exclude '.' & '..'
    
        mouseName = fileStruct(f).name(1:5);
        targetPath = ['C:\Users\Mate Neubrandt\UIO Physiology Dropbox Dropbox\Mate Neubrandt\MateData\RECORDINGS' '\' mouseName];
        if ~isfolder(targetPath); mkdir(targetPath); end
        sourceFile = fullfile(fileStruct(f).folder,fileStruct(f).name);
        destinationFile = fullfile(targetPath,fileStruct(f).name);
        
        copyfile(sourceFile,destinationFile)
        delete(sourceFile)

        targetPaths{f} = targetPath;

    end
end
disp('Files have been copied to the target folder.')


