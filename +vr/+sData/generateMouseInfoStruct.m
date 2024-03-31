%% Generate or edit mouseInfo structure 
%
% Keep all mousinfo file in a spectfied folder, use the mouse name as
% filename. These mouse info files will be imported when generating sData
% files.
% Keep a template struct with the file name "default" in the same
% directory.

clear;

% Set paramaters
mouseInfoFolder = 'C:\Users\Mate Neubrandt\UIO Physiology Dropbox Dropbox\Mate Neubrandt\MateData\RECORDINGS\MOUSEINFO';


%% Select existing mousinfo file 

[fileName,~,~] = uigetfile('*.mat','Select existing file to edit or as a template. Press "Cancel" to start with default template.',mouseInfoFolder,'MultiSelect','off' );

if fileName == 0
    fileName = 'default';
end

load(fullfile(mouseInfoFolder,fileName));


% Extract values from default struct 
name = mouseInfo.name;
dateOfBirth = mouseInfo.dateOfBirth;
strain = mouseInfo.strain;
transgene = mouseInfo.transgene;
link = mouseInfo.linkToVendor;
sex = mouseInfo.sex;
lightCycle = mouseInfo.lightCycle;
SLCageNumber = num2str(mouseInfo.SLCageNumber);
SLMouseNumber = num2str(mouseInfo.SLMouseNumber);
surgeryDate = mouseInfo.surgeryDate;
windowCoordinates = num2str(mouseInfo.windowCoordinates);
surgeryDoneBy = mouseInfo.surgeryDoneBy;
surgeryProtocol = mouseInfo.surgeryProtocol;
windowType = mouseInfo.windowType;
injectedVirus = mouseInfo.injectedVirus;
injectedVirusLocation = mouseInfo.injectedVirusLocation;
if isa(mouseInfo.injectedVirusNanoLPerSite,'double')
    injectedVirusNanoLPerSite = num2str(mouseInfo.injectedVirusNanoLPerSite);
elseif ischar(mouseInfo.injectedVirusNanoLPerSite)
    injectedVirusNanoLPerSite = mouseInfo.injectedVirusNanoLPerSite;
end
if isa(mouseInfo.injectedVirusNumberOfSites,'double')
    injectedVirusNumberOfSites = num2str(mouseInfo.injectedVirusNumberOfSites);
elseif ischar(mouseInfo.injectedVirusNumberOfSites)
    injectedVirusNumberOfSites = mouseInfo.injectedVirusNumberOfSites;
end
if isa(mouseInfo.injectedVirusDepth,'double')
    injectedVirusDepth = num2str(mouseInfo.injectedVirusDepth);
elseif ischar(mouseInfo.injectedVirusDepth)
    injectedVirusDepth = mouseInfo.injectedVirusDepth;
end




% Manual input of some data
prompt = {'Mouse name:','Date of birth (YYY.MM.DD):','Strain:','Transgene:','Webpage link for transgene:',...
    'Sex','Light cycle:','Science Linker cage number:','Science Linker mouse number:'};

dlgtitle = 'Fill In These Data';
definput = {name,dateOfBirth,strain,transgene,link,...
    sex,lightCycle,SLCageNumber,SLMouseNumber};
answer = inputdlg(prompt,dlgtitle,[1 40],definput);



% Update values 
name = answer{1};
dateOfBirth = answer{2};
strain = answer{3};
transgene = answer{4};
link = answer{5};
sex = answer{6};
lightCycle = answer{7};
SLCageNumber = answer{8};
SLMouseNumber = answer{9};



% Manual input of some data
prompt = {'Surgery date:','Window coordinates (X Y, separate by space):','Who did the surgery?','Surgery protocol:','Window type:',...
    'Injected virus:','Location of virus injection:','Amount of injected virus (nl/site)','Number of injection sites:','Depth of injection:'};

dlgtitle = 'Fill In These Data';
definput = {surgeryDate,windowCoordinates,surgeryDoneBy,surgeryProtocol,windowType,...
    injectedVirus,injectedVirusLocation,injectedVirusNanoLPerSite,injectedVirusNumberOfSites,injectedVirusDepth};
opts.Resize = 'on';
answer2 = inputdlg(prompt,dlgtitle,[1 40],definput);



% Update values 
surgeryDate = answer2{1};
windowCoordinates = answer2{2};
surgeryDoneBy = answer2{3};
surgeryProtocol = answer2{4};
windowType = answer2{5};
injectedVirus = answer2{6};
injectedVirusLocation = answer2{7};
injectedVirusNanoLPerSite = answer2{8};
injectedVirusNumberOfSites = answer2{9};
injectedVirusDepth = answer2{10};


%% Save values to struct 
mouseInfo.name = name;
mouseInfo.dateOfBirth = dateOfBirth;
mouseInfo.strain = strain;
mouseInfo.transgene = transgene;
mouseInfo.linkToVendor = link;
mouseInfo.sex = sex;
mouseInfo.lightCycle = lightCycle;

mouseInfo.SLCageNumber = str2double(SLCageNumber);
mouseInfo.SLMouseNumber = str2double(SLMouseNumber);
mouseInfo.surgeryDate = surgeryDate;
mouseInfo.windowCoordinates = str2num(windowCoordinates);
mouseInfo.surgeryDoneBy = surgeryDoneBy;
mouseInfo.surgeryProtocol = surgeryProtocol;
mouseInfo.windowType = windowType;
mouseInfo.injectedVirus = injectedVirus;
mouseInfo.injectedVirusLocation = injectedVirusLocation;

if numel(str2num(injectedVirusNanoLPerSite)) > 0
    mouseInfo.injectedVirusNanoLPerSite = str2num(injectedVirusNanoLPerSite);
else
    mouseInfo.injectedVirusNanoLPerSite = injectedVirusNanoLPerSite;
end

if numel(str2num(injectedVirusNumberOfSites)) > 0
    mouseInfo.injectedVirusNumberOfSites = str2num(injectedVirusNumberOfSites);
else
    mouseInfo.injectedVirusNumberOfSites = injectedVirusNumberOfSites;
end

if numel(str2num(injectedVirusDepth)) > 0
    mouseInfo.injectedVirusDepth = str2num(injectedVirusDepth);
else
    mouseInfo.injectedVirusDepth = injectedVirusDepth;
end


%% Save new file or overwrite existing file in case of editing

save(fullfile(mouseInfoFolder,mouseInfo.name),'mouseInfo');





%%

%{
% Manual input of some data
prompt = {'Mouse name:','Date of birth:','Strain:','Transgene:','Webpage link for transgene:',...
    'Sex','Light cycle:','Science Linker cage number:','Science Linker mouse number:','Surgery date:',...
    'Window coordinate X:','Window coordinate Y:','Who did the surgery?','Surgery protocol:','Window type:',...
    'Injected virus:','Location of virus injection:','Amount of injected virus (nl/site)','Number of injection sites:','Depth of injection:'};

dlgtitle = 'Fill In These Data';

definput = {name,dateOfBirth,strain,transgene,ordered,...
    sex,lightCycle,SLCageNumber,SLMouseNumber,surgeryDate,...
    windowCoordinatesX,windowCoordinatesY,surgeryDoneBy,surgeryProtocol,windowType,...
    injectedVirus,injectedVirusLocation,injectedVirusNanoLPerSite,injectedVirusNumberOfSites,injectedVirusDepth};
opts.Resize = 'on';
answer = inputdlg(prompt,dlgtitle,[1 40],definput);


%}