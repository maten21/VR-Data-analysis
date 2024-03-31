% Plot opto figure median and quartiles
% [fileName filePath] = uigetfile('*.mat');
function [] = opto2Figs(fileName,filePath)

load(fullfile(filePath,fileName));


%% Settings

nTrialTypes = 2; %Use it to derermine which stimulus prorocol to plot in addition to the control
smoothSpan = 5;
plotLegend = {'no stim'; ''; 'Full trial'; ''};
rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;

mymap = lines(nTrialTypes);
Xax = sData.stats.sessionAvs(1).plotXAxis;
binVel = sData.behavior.trialMatrices.binVel;
lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;

lickDistrCtrl = sData.trials.trialTypesMeta(1).lickQuartiles(2)*rewardZone/100;
lickDistrStim = sData.trials.trialTypesMeta(nTrialTypes).lickQuartiles(2)*rewardZone/100;
firstLickCtrl = sData.trials.trialTypesMeta(1).firstLicksCm(2);
firstLickStim = sData.trials.trialTypesMeta(nTrialTypes).firstLicksCm(2);

% for plotting stimulus location
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(nTrialTypes).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));
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
rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binVel(sData.trials.trialTypesMeta(t).trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
end


axis([corridorStart corridorLength 0 70]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.0025,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    
    color = mymap(t,:);
    trialMatrix = lickFreq(sData.trials.trialTypesMeta(t).trials,:);
    
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
end


rectangle('Position',[lickDistrCtrl-1, 4.0025, 2, 10],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[lickDistrStim-1, 4.0025, 2, 10],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

rectangle('Position',[firstLickCtrl-0.5, 0.0025, 1, 4],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[firstLickStim-0.5, 0.0025, 1, 4],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
legend(plotLegend,'Location','northwest')



if ~exist([filePath 'optoFigs'], 'dir')
    mkdir([filePath 'optoFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'optoFigs',[sData.sessionInfo.sessionID '_Quartiles' '-trialType' num2str(nTrialTypes)]),'.png'));
close(gcf)




%% Settings

nTrialTypes = 3; %Use it to derermine which stimulus prorocol to plot in addition to the control
smoothSpan = 5;
plotLegend = {'no stim'; ''; 'Home-box'; ''};
rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;

mymap = lines(nTrialTypes);
Xax = sData.stats.sessionAvs(1).plotXAxis;
binVel = sData.behavior.trialMatrices.binVel;
lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;

lickDistrCtrl = sData.trials.trialTypesMeta(1).lickQuartiles(2)*rewardZone/100;
lickDistrStim = sData.trials.trialTypesMeta(nTrialTypes).lickQuartiles(2)*rewardZone/100;
firstLickCtrl = sData.trials.trialTypesMeta(1).firstLicksCm(2);
firstLickStim = sData.trials.trialTypesMeta(nTrialTypes).firstLicksCm(2);

% for plotting stimulus location
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(nTrialTypes).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));
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
rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binVel(sData.trials.trialTypesMeta(t).trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
end


axis([corridorStart corridorLength 0 70]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.0025,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    
    color = mymap(t,:);
    trialMatrix = lickFreq(sData.trials.trialTypesMeta(t).trials,:);
    
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
end


rectangle('Position',[lickDistrCtrl-1, 4.0025, 2, 10],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[lickDistrStim-1, 4.0025, 2, 10],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

rectangle('Position',[firstLickCtrl-0.5, 0.0025, 1, 4],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[firstLickStim-0.5, 0.0025, 1, 4],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
legend(plotLegend,'Location','northwest')



if ~exist([filePath 'optoFigs'], 'dir')
    mkdir([filePath 'optoFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'optoFigs',[sData.sessionInfo.sessionID '_Quartiles' '-trialType' num2str(nTrialTypes)]),'.png'));
close(gcf)


%% Settings

nTrialTypes = 4; %Use it to derermine which stimulus prorocol to plot in addition to the control
smoothSpan = 5;
plotLegend = {'no stim'; ''; 'On-the-path'; ''};
rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;

mymap = lines(nTrialTypes);
Xax = sData.stats.sessionAvs(1).plotXAxis;
binVel = sData.behavior.trialMatrices.binVel;
lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;

lickDistrCtrl = sData.trials.trialTypesMeta(1).lickQuartiles(2)*rewardZone/100;
lickDistrStim = sData.trials.trialTypesMeta(nTrialTypes).lickQuartiles(2)*rewardZone/100;
firstLickCtrl = sData.trials.trialTypesMeta(1).firstLicksCm(2);
firstLickStim = sData.trials.trialTypesMeta(nTrialTypes).firstLicksCm(2);

% for plotting stimulus location
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(nTrialTypes).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));
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
rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binVel(sData.trials.trialTypesMeta(t).trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
end


axis([corridorStart corridorLength 0 70]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.0025,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    
    color = mymap(t,:);
    trialMatrix = lickFreq(sData.trials.trialTypesMeta(t).trials,:);
    
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
end


rectangle('Position',[lickDistrCtrl-1, 4.0025, 2, 10],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[lickDistrStim-1, 4.0025, 2, 10],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

rectangle('Position',[firstLickCtrl-0.5, 0.0025, 1, 4],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[firstLickStim-0.5, 0.0025, 1, 4],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
legend(plotLegend,'Location','northwest')



if ~exist([filePath 'optoFigs'], 'dir')
    mkdir([filePath 'optoFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'optoFigs',[sData.sessionInfo.sessionID '_Quartiles' '-trialType' num2str(nTrialTypes)]),'.png'));
close(gcf)

%{

nTrialTypes = 5; %Use it to derermine which stimulus prorocol to plot in addition to the control
smoothSpan = 5;
plotLegend = {'no stim'; ''; '10 - reward'; ''};
rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;

mymap = lines(nTrialTypes);
Xax = sData.stats.sessionAvs(1).plotXAxis;
binVel = sData.behavior.trialMatrices.binVel;
lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;

lickDistrCtrl = sData.trials.trialTypesMeta(1).lickQuartiles(2)*rewardZone/100;
lickDistrStim = sData.trials.trialTypesMeta(nTrialTypes).lickQuartiles(2)*rewardZone/100;
firstLickCtrl = sData.trials.trialTypesMeta(1).firstLicksCm(2);
firstLickStim = sData.trials.trialTypesMeta(nTrialTypes).firstLicksCm(2);

% for plotting stimulus location
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(nTrialTypes).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));
% Generate plot legend cell array AUTOMATIC
%{
for t = 1:1:nTrialTypes
    vr.trialType;
    protType = strsplit(sData.trials.trialTypesMeta(t).stimProtocol,' int');

    plotLegend{t} = [protType{1} ' (n = ' num2str(nTrials) ')'];
end
%}




figure('Color','white','Position',[0 0 850 350]); %pos of figure [left bottom width height]

suptitle([sData.sessionInfo.sessionID newline])

subplot(1,2,1); % Velocity
hold on

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binVel(sData.trials.trialTypesMeta(t).trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
end


axis([corridorStart corridorLength 0 70]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.0025,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    
    color = mymap(t,:);
    trialMatrix = lickFreq(sData.trials.trialTypesMeta(t).trials,:);
    
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
end


rectangle('Position',[lickDistrCtrl-1, 4.0025, 2, 10],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[lickDistrStim-1, 4.0025, 2, 10],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

rectangle('Position',[firstLickCtrl-0.5, 0.0025, 1, 4],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[firstLickStim-0.5, 0.0025, 1, 4],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
legend(plotLegend,'Location','northwest')



if ~exist([filePath 'optoFigs'], 'dir')
    mkdir([filePath 'optoFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'optoFigs',[sData.sessionInfo.sessionID '_Quartiles' '-trialType' num2str(nTrialTypes)]),'.png'));
close(gcf)


%}


%% Mean std

%{
figure('Color','white','Position',[0 0 850 350]); %pos of figure [left bottom width height]

suptitle([sData.sessionInfo.sessionID newline])

subplot(1,2,1); % Velocity
hold on

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binVel(sData.trials.trialTypesMeta(t).trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotMeanStd(Xax,trialMatrix,color,smoothSpan)
end

axis([corridorStart corridorLength 0 70]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.0025,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus

%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

map = lines(nTrialTypes);
for t = 1:nTrialTypes-1:nTrialTypes
    
    color = mymap(t,:);
    trialMatrix = lickFreq(sData.trials.trialTypesMeta(t).trials,:);
    
    plotMeanStd(Xax,trialMatrix,color,smoothSpan)
    
end

%rectangle('Position',[lickDistrCtrl-1, 0.0025, 2, 80],'FaceColor',mymap(1,:),'EdgeColor','none'); % home-box
%rectangle('Position',[lickDistrStim-1, 0.0025, 2, 80],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % RZ

axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
legend(plotLegend,'Location','northwest')


saveas(gcf,strcat(fullfile(filePath,'optoFigs',[sData.sessionInfo.sessionID '_2-meanStd' 'trialType' num2str(nTrialTypes)]),'.png'));
close(gcf)
%}

end