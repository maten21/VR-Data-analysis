% Plot opto figure median and quartiles
% [fileName filePath] = uigetfile('*.mat');

function [] = opto2Figs(fileName,filePath)

if nargin == 2
    load(fullfile(filePath,fileName));
else
    [fileName filePath] = uigetfile('*.mat');
    load(fullfile(filePath,fileName));
end
    

    
%% Settings

%{
switch sData.trials.trialTypesMeta(i).trialTypeIndicator
    case 1
        stimulusType = 2;
        plotLegend = {'no stim'; ''; 'Full trial'; ''};
    case 2
        stimulusType = 3;    
        plotLegend = {'no stim'; ''; 'Home-box'; ''};
    case 5
        stimulusType = 4;
        plotLegend = {'no stim'; ''; 'On-the-path'; ''};
    otherwise
        msgbox('ERROR: Ttial type indicator is not configured!')
end
%}

nTrialTypes = size(sData.trials.trialTypesMeta,2);
lineWidth = 1.5;
plotLegend = {'no stim'; 'Home-box'; 'On-the-path'};
rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;

mymap = lines(nTrialTypes);
Xax = sData.stats.sessionAvs(1).plotXAxis;
binVel = sData.behavior.trialMatrices.binVel;
lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;

%{
lickDistrCtrl = sData.trials.trialTypesMeta(1).lickQuartiles(2)*rewardZone/100;
lickDistrStim = sData.trials.trialTypesMeta(i).lickQuartiles(2)*rewardZone/100;
firstLickCtrl = sData.trials.trialTypesMeta(1).firstLicksCm(2);
firstLickStim = sData.trials.trialTypesMeta(i).firstLicksCm(2);

% for plotting stimulus location
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(i).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));
%}

% Generate plot legend cell array AUTOMATIC
%{
for t = 1:1:nTrialTypes
    vr.trialType;
    protType = strsplit(sData.trials.trialTypesMeta(t).stimProtocol,' int');

    plotLegend{t} = [protType{1} ' (n = ' num2str(nTrials) ')'];
end
%}

%% Plot median figure

figure('Color','white','Position',[0 0 850 350]); %pos of figure [left bottom width height]

suptitle([sData.sessionInfo.sessionID newline])

subplot(1,2,1); % Velocity
hold on

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ



for i = 1:1:nTrialTypes   
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(i).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));

color = mymap(i,:);
if i>1
    
    rectangle('Position',[stimStart,76-((i-2)*4),(stimEnd-stimStart),4],'FaceColor',color,'EdgeColor','none'); % Stimulus
end

%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone


 
    %color = mymap(1,:);
    %trialMatrix = binVel(sData.trials.trialTypesMeta(1).trials,:);
    %plotQuartiles(Xax,trialMatrix,color,smoothSpan)

    %color = mymap(i,:);
    trialMatrix = binVel(sData.trials.trialTypesMeta(i).trials,:);
    plot(Xax,nanmedian(trialMatrix),'Color',color,'LineWidth',lineWidth);
    
end


axis([corridorStart corridorLength 0 80]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.05,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ

for i = 1:1:nTrialTypes   
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(i).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));

color = mymap(i,:);
if i>1    
    rectangle('Position',[stimStart,7.6-((i-2)*0.4),(stimEnd-stimStart),0.4],'FaceColor',color,'EdgeColor','none'); % Stimulus
end
firstLick = sData.trials.trialTypesMeta(i).firstLicksCm(2);
rectangle('Position',[firstLick-1, 6.0025, 2, 4],'FaceColor',color,'EdgeColor','none'); % stim lick distr


%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
    
    %color = mymap(1,:);
    %trialMatrix = lickFreq(sData.trials.trialTypesMeta(1).trials,:); 
    %plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
    %color = mymap(i,:);
    trialMatrix = lickFreq(sData.trials.trialTypesMeta(i).trials,:); 
    plot(Xax,nanmedian(trialMatrix),'Color',color,'LineWidth',lineWidth);
    
end    


axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
legend(plotLegend,'Location','northwest')




saveas(gcf,strcat(fullfile(filePath,[sData.sessionInfo.sessionID '_optoFig']),'.png'));
close(gcf)




end