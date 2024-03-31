



sessionFolderRawdata = fullfile(sessionObject.DataLocation(1).RootPath,sessionObject.DataLocation(1).Subfolders);
rootFolder = fileparts(sessionFolderRawdata); rootFolder = fileparts(rootFolder); rootFolder = fileparts(rootFolder);
expDateFolder = fileparts(sessionFolderRawdata);
sessionFolderProcessed = sessionObject.getSessionFolder();



optionsText = {'Create cranial window', 'Align scanfields', 'Plot FOVs'};

if isfile(fullfile(rootFolder,'windowInventory',[sessionObject.subjectID '.mat']))
    load(fullfile(rootFolder,'windowInventory',[sessionObject.subjectID '.mat']))
    if length(fovDb.Windows) > 0

        if isfile(fullfile(sessionFolderProcessed,'FovManagerData',[sessionObject.sessionID '_scanArea.mat']))
            optionsText = optionsText(1:3);
            initVal = 3;

            %if isfile(fullfile(sessionFolderProcessed,'FovManagerData',[sessionObject.sessionID '_FOVs.mat']))
            %end

        else % Align scan area first
            optionsText = optionsText(1:2);
            initVal = 2;
        end
    else % No window is saved. Create window first.
        optionsText = optionsText(1);
        initVal = 1;
    end

else % Create window first
    optionsText = optionsText(1);
    initVal = 1;
end


% Dialoge box to select what to do
selection = listdlg('PromptString',{'FOV data can be edited on', 'different hierarchical levels.', 'Select what to do and save', 'inventory when finished.'},...
    'SelectionMode','single','ListString',optionsText,'InitialValue',initVal);




%% Create window

if selection == 1

if isfile(fullfile(rootFolder,'windowInventory',[sessionObject.subjectID '.mat']))
    load(fullfile(rootFolder,'windowInventory',[sessionObject.subjectID '.mat']))
else
    load(fullfile(rootFolder,'windowInventory',['inventoryTemplate', '.mat']))
    fovDb.MouseId = sessionObject.subjectID;
    loadPath = fullfile(rootFolder,'windowInventory',[sessionObject.subjectID, '.mat']);
    save(loadPath,"fovDb");
end



hFovmanager = fovmanager;
hFovmanager.loadFovDatabase(loadPath)

end



%% Align scanfields

if selection == 2

%if isfile(fullfile(rootFolder,'windowInventory',[sessionObject.subjectID '.mat']))
    load(fullfile(rootFolder,'windowInventory',[sessionObject.subjectID '.mat']))
%end

[fullFovFileNames, fullFovFilePaths] = vr.findFileInFolder(expDateFolder, sessionObject.subjectID, 'tif');


for i = 1:1:length(fullFovFilePaths)

    
    im = vr.mesoscope.readMultiFovTiff(fullFovFilePaths{i});
    fullFovData = vr.mesoscope.ScanImageFovData(fullFovFilePaths{i});

    im = mean(im,3);
    lutMin = quantile(im(:),0.01);
    lutMax = quantile(im(:),0.99);

    im = (im-lutMin)./lutMax;
    im = im2uint8(im); % convert to 8bit

    im = imresize(im, 0.5);



% Calculate scan field file size

for f = 1:1:length(fullFovData)

    centerX(f) = fullFovData(f).fovCenterUmXY(1);
    centerY(f) = fullFovData(f).fovCenterUmXY(2);
    sizeX(f) = fullFovData(f).fovSizeUmXY(1);
    sizeY(f) = fullFovData(f).fovSizeUmXY(2);


end

[maxCenterX, iMaxCenterX] = max(centerX);
[minCenterX, iMinCenterX] = min(centerX);
[maxCenterY, iMaxCenterY] = max(centerY);
[minCenterY, iMinCenterY] = min(centerY);

figSizeX = maxCenterX - minCenterX + sizeX(iMinCenterX)/2 + sizeX(iMaxCenterX)/2;
figSizeY = maxCenterY - minCenterY + sizeY(iMinCenterY)/2 + sizeY(iMaxCenterY)/2;

Ydir = -1; % 1 or -1

corners = [-figSizeX/2, -figSizeY/2*Ydir; ...
           -figSizeX/2, figSizeY/2*Ydir; ...
            figSizeX/2, figSizeY/2*Ydir; ...
            figSizeX/2, -figSizeY/2*Ydir];

corners = corners/1000;



% Generate fovArray

%fovArray = struct();
fovArray(i).alpha = 0.3; 
fovArray(i).center = [0, 0]; 
fovArray(i).currentSession = [];
fovArray(i).depth = [];
fovArray(i).edge = corners;
fovArray(i).image = im;
fovArray(i).listOfSessions = [];
fovArray(i).nRois = [];
fovArray(i).orientation.isMirroredX = false;
fovArray(i).orientation.isMirroredY = true;
fovArray(i).orientation.theta = [];
fovArray(i).shape = 'square';
% fovArray(i).imagingSystem = fullFovData(1).imagingSystem; %index is intentionally 1


end

fovDb.MouseId = sessionObject.subjectID;
fovDb.Windows.alpha = 0.7;
fovDb.Windows.fovArray = fovArray;


if ~isfolder(fullfile(sessionFolderProcessed,'FovManagerData')); mkdir(fullfile(sessionFolderProcessed,'FovManagerData')); end

loadPath = fullfile(sessionFolderProcessed,'FovManagerData',[sessionObject.sessionID '_scanArea.mat']);
save(loadPath,"fovDb");

hFovmanager = fovmanager;
hFovmanager.loadFovDatabase(loadPath)



end


%% Load exp FOVs 

if selection == 3


load(fullfile(rootFolder,'windowInventory',[sessionObject.subjectID '.mat']))

fovDbScanField = load(fullfile(sessionFolderProcessed,'FovManagerData',[sessionObject.sessionID '_scanArea.mat']));



fileName = [sprintf('%s_XYT', sessionObject.sessionID), '_meta2P', '.mat'];
filePath = fullfile(sessionFolderProcessed,'metadata');


if isfile(fullfile(filePath,fileName))
    load(fullfile(filePath,fileName));
else
    msgbox(['Metadata file is not found: ' fullfile(filePath,fileName)])
end

fullFovFileNames = vr.findFileInFolder(expDateFolder, sessionObject.subjectID, 'tif'); % check name for path ID

[~, alignedDataFolders] = vr.findFileInFolder(sessionFolderProcessed, {'FOV', 'fov'}, 'dir');


for i = 1:1:length(meta2P.fovData)



    for g = 1:1:length(alignedDataFolders)
        if isfolder(alignedDataFolders{g}) && str2num(alignedDataFolders{g}(end)) == i
            path = fullfile(alignedDataFolders{g},'\fov_images');
            files = vr.findFileInFolder(path, 'average_projection_corr.tif');
            im = tiffreadVolume(fullfile(path, files{1}));
            % hImFov = imshow(fovIm); %,[lutMin, lutMax]);

        end
    end

figSizeX = meta2P.fovData(i).fovSizeUmXY(1);
figSizeY = meta2P.fovData(i).fovSizeUmXY(2);

centerX = meta2P.fovData(i).fovCenterUmXY(1);
centerY = meta2P.fovData(i).fovCenterUmXY(2);


for f = 1:1:length(fullFovFileNames)

    if numel(strfind(fullFovFileNames{f},'_P1_')) > 0
        origoPos = fovDbScanField.fovDb.Windows.fovArray(1).center*1000;
    elseif numel(strfind(fullFovFileNames{f},'_P2_')) > 0
        origoPos = fovDbScanField.fovDb.Windows.fovArray(2).center*1000;
    end
end

Ydir = -1; % 1 or -1

centerX = centerX + origoPos(1);
centerY = centerY*Ydir + origoPos(2);



corners = [centerX - figSizeX/2, centerY - figSizeY/2*Ydir; ...
           centerX - figSizeX/2, centerY + figSizeY/2*Ydir; ...
           centerX + figSizeX/2, centerY + figSizeY/2*Ydir; ...
           centerX + figSizeX/2, centerY - figSizeY/2*Ydir];

corners = corners/1000;




%% Generate fovArray

%fovArray = struct();

fovArray(i).alpha = 1; 
fovArray(i).center = [0, 0]; %[centerX, centerY]/1000; 
fovArray(i).currentSession = [];
fovArray(i).depth = [];
fovArray(i).edge = corners;
fovArray(i).image = im;
fovArray(i).listOfSessions = [];
fovArray(i).nRois = [];
fovArray(i).orientation.isMirroredX = false;
fovArray(i).orientation.isMirroredY = true;
fovArray(i).orientation.theta = [];
fovArray(i).shape = 'square';



end

%fovDb.MouseId = sessionObject.subjectID;
fovDb.Windows.alpha = 0.7;
fovDb.Windows.fovArray = fovArray;


if ~isfolder(fullfile(sessionFolderProcessed,'FovManagerData')); mkdir(fullfile(sessionFolderProcessed,'FovManagerData')); end

loadPath = fullfile(sessionFolderProcessed,'FovManagerData',[sessionObject.sessionID '_FOVs.mat']);

save(loadPath,"fovDb");

hFovmanager = fovmanager;
hFovmanager.loadFovDatabase(loadPath)


end


