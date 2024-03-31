clear;

[sDataFiles,filePath] = vr.loadData('light');

nFiles = size(sDataFiles,2);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correct trial indicator errors!!!

%% Calculate summary curves

if exist('curves')
    clear('curves')
end
curves(15).vel = 0;
curves(15).freq = 0;
curves(15).optStim = 0;
curves(15).prot = '';

for f = 1:1:nFiles
    
    for i = 1:1:size(sDataFiles{1,f}.trials.trialTypesMeta,2)
        t = sDataFiles{1,f}.trials.trialTypesMeta(i).trialTypeIndicator + 1;
        r = size(curves(t).freq,1) + 1;
        
        curves(t).vel(r,:) = sDataFiles{1,f}.stats.sessionMedians(i).medBinVel;
        curves(t).freq(r,:) = sDataFiles{1,f}.stats.sessionMedians(i).medLickFreqInBin;
        curves(t).optStim(r,:) = sDataFiles{1,f}.stats.sessionAvs(i).avOptStimMatrix;
        curves(t).prot{r,:} = sDataFiles{1,f}.trials.trialTypesMeta(i).stimProtocol;
        
    end
    
end

%changed to median this is the original with mean
%{ 
for t = 1:1:10
    if size(curves(t).vel,1)>0
        avCurves(t).vel = mean(curves(t).vel);
        avCurves(t).freq = mean(curves(t).freq);
        avCurves(t).optStim = mean(curves(t).optStim);
        avCurves(t).prot = curves(t).prot{1,1};
    end
end
%}

for t = 1:1:10
    if size(curves(t).vel,1)>0
        avCurves(t).vel = median(curves(t).vel);
        avCurves(t).freq = median(curves(t).freq);
        avCurves(t).optStim = mean(curves(t).optStim);
        avCurves(t).prot = curves(t).prot{1,1};
    end
end

% for plotting stimulus location
for i = 1:1:numel(avCurves) % numel(subset)
    t = i;
    % t = subset(i);
    %stimCurves(t) = nanmean(optStimMatrix(sData.trials.trialTypesMeta(nTrialTypes).trials,:));
    stimLocation(t).optStim = avCurves(t).optStim;
    
    stimMax = max(avCurves(t).optStim);
    stimMin = min(avCurves(t).optStim);
    if stimMax < stimMin*1.1 % control
        stimLocation(t).optStim = zeros(size(stimLocation(t).optStim));
    else
        stimLocation(t).optStim(avCurves(t).optStim < (stimMax - stimMin)/3) = 0;
        stimLocation(t).optStim(avCurves(t).optStim >= (stimMax - stimMin)/3) = 1;
    end
end


% PLOT CURVES
subset = [1, 2, 3, 6]; % full trial, homebox, OTP
subset = [1, 3, 6]; %  homebox, OTP
subset = [2, 3, 4, 5, 6, 7, 8, 9, 10];
subset = [2, 3, 4, 5, 6, 7, 10];
subset = [2, 3, 8, 9];
subset = [2, 3, 4, 6];
subset = [1, 2, 3, 6];
subset = [3, 4, 6]; % home box, 60-120, 10-70
subset = [5, 6, 7]; % diff lengths
subset = [4, 6, 7, 10]; % OTP
subset = [4, 5, 6, 7, 10]; % OTP all

rewardZone = sDataFiles{1, 1}.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;
viewDistance = sDataFiles{1, 1}.behavior.viewDistance;  

mymap = [lines(7); 0.5 0.5 0.5; 1 0.1 0.1; 0 1 0.5 ];

plotXAxis = sDataFiles{1, 1}.stats.sessionAvs(1).plotXAxis;




figure('Color','white','Position',[0 0 850 350]); %pos of figure [left bottom width height]

lineWidth = 1.5;
%suptitle([sData.sessionInfo.sessionID newline])

subplot(1,2,1); % Velocity
hold on
adjustY = 15;

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80+adjustY],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80+adjustY],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
rectangle('Position',[-1-viewDistance,0.5,2,80+adjustY],'FaceColor',[0.9 0.9 0.97],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i = 1:1:numel(subset)
    t = subset(i);
    %vr.trialType;
    color = mymap(t,:);
    data = avCurves(t).vel;
    plot(plotXAxis,data','Color',color,'LineWidth',lineWidth)
    
    if t > 1
        data = stimLocation(t).optStim(stimLocation(t).optStim > 0)*((70+adjustY)-(70+adjustY)/50*t);
        plot(plotXAxis(stimLocation(t).optStim > 0),data,'Color',color,'LineWidth',lineWidth*2)
    end
end
text(0,70+adjustY,'\fontsize{11}no stim','Color',mymap(1,:))

axis([corridorStart corridorLength 0 70+adjustY]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.05,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.05,10,80+adjustY],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
rectangle('Position',[-1-viewDistance,0.05,2,80],'FaceColor',[0.9 0.9 0.97],'EdgeColor','none'); % RZ
%rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i = 1:1:numel(subset)
    t = subset(i);
    
    color = mymap(t,:);
    data = avCurves(t).freq;
    plot(plotXAxis,data,'Color',color,'LineWidth',lineWidth)
    
    if t > 1
        data = stimLocation(t).optStim(stimLocation(t).optStim > 0)*(8-8/50*t);
        plot(plotXAxis(stimLocation(t).optStim > 0),data,'Color',color,'LineWidth',lineWidth*2)
    end
end
text(0,8,'\fontsize{11}no stim','Color',mymap(1,:))


axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick rate (Hz)');
%legend(plotLegend,'Location','northwest')




%% SAVE FIG
% fileName = 'ctrlStandard'; 
% fileName = 'ctrlDiffLengths'; 
% fileName = 'controlSummary1'; 
% fileName = 'controlSummary2'; 
% fileName = 'opsinStandard'; 
% fileName = 'opsinStandard6mice2SessionsNorm'; 
% fileName = 'opsinStandard6Mice3SessionsMed';
% fileName = 'opsinOTPs'; 
% fileName = 'opsinDiffLengths'; 
% fileName = 'opsinShortStim'; 
% fileName = 'opsinSummary'; 
% fileName = 'opsinSummary2'; 
% fileName = 'maskStimSummary'; 

if ~exist([filePath 'summFigs'], 'dir')
    mkdir([filePath 'summFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'summFigs',fileName),'.png'));
close(gcf)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate NORMALIZED summary curves

if exist('curves')
    clear('curves')
end
curves(15).vel = 0;
curves(15).freq = 0;
curves(15).optStim = 0;
curves(15).prot = '';

for f = 1:1:nFiles
    
    for i = 2:1:size(sDataFiles{1,f}.trials.trialTypesMeta,2)
        t = sDataFiles{1,f}.trials.trialTypesMeta(i).trialTypeIndicator + 1;
        r = size(curves(t).freq,1) + 1;
        
        curves(t).vel(r,:) = sDataFiles{1,f}.stats.sessionMedians(i).medBinVel - sDataFiles{1,f}.stats.sessionMedians(1).medBinVel;
        curves(t).freq(r,:) = sDataFiles{1,f}.stats.sessionMedians(i).medLickFreqInBin - sDataFiles{1,f}.stats.sessionMedians(1).medLickFreqInBin;
        curves(t).optStim(r,:) = sDataFiles{1,f}.stats.sessionAvs(i).avOptStimMatrix;
        curves(t).prot{r,:} = sDataFiles{1,f}.trials.trialTypesMeta(i).stimProtocol;
        
    end
    
end

for t = 1:1:10
    if size(curves(t).vel,1)>0
        avCurves(t).vel = mean(curves(t).vel);
        avCurves(t).freq = mean(curves(t).freq);
        avCurves(t).optStim = mean(curves(t).optStim);
        avCurves(t).prot = curves(t).prot{1,1};
    end
end

%{
figure
hold on
plotXAxis = sDataFiles{1, 1}.stats.sessionAvs(1).plotXAxis;
plot(plotXAxis,avCurves(1).freq)
plot(plotXAxis,avCurves(2).freq)
plot(plotXAxis,avCurves(3).freq)
plot(plotXAxis,avCurves(4).freq)

figure
hold on
plotXAxis = sDataFiles{1, 1}.stats.sessionAvs(1).plotXAxis;
plot(plotXAxis,avCurves(1).vel)
plot(plotXAxis,avCurves(2).vel)
plot(plotXAxis,avCurves(3).vel)
plot(plotXAxis,avCurves(4).vel)

figure
hold on
plotXAxis = sDataFiles{1, 1}.stats.sessionAvs(1).plotXAxis;
plot(plotXAxis,avCurves(1).optStim)
plot(plotXAxis,avCurves(2).optStim)
plot(plotXAxis,avCurves(3).optStim)
plot(plotXAxis,avCurves(4).optStim)

%}


%{
%%%%%%%% Previous working version

rewardZone = sDataFiles{1, 1}.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;
viewDistance = sDataFiles{1, 1}.behavior.viewDistance;  

mymap = lines(nTrialTypes);
plotXAxis = sDataFiles{1, 1}.stats.sessionAvs(1).plotXAxis;
%binVel = sData.behavior.trialMatrices.binVel;
%lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
%optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;

% for plotting stimulus location
for t = 1:1:nTrialTypes
    %stimCurves(t) = nanmean(optStimMatrix(sData.trials.trialTypesMeta(nTrialTypes).trials,:));
    stimLocation(t).optStim = avCurves(t).optStim;
    
    stimMax = max(avCurves(t).optStim);
    stimMin = min(avCurves(t).optStim);
    if stimMax < stimMin*1.1 % control
        stimLocation(t).optStim = zeros(size(stimLocation(t).optStim));
    else
        stimLocation(t).optStim(avCurves(t).optStim < (stimMax - stimMin)/3) = 0;
        stimLocation(t).optStim(avCurves(t).optStim >= (stimMax - stimMin)/3) = 1;
    end
end

% stimStart = plotXAxis(min(find(stimLocation == 1)));
% stimEnd = Xax(max(find(stimLocation == 1)));


figure('Color','white','Position',[0 0 850 350]); %pos of figure [left bottom width height]

lineWidth = 1.5;
%suptitle([sData.sessionInfo.sessionID newline])

subplot(1,2,1); % Velocity
hold on

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
rectangle('Position',[-1-viewDistance,0.5,2,80],'FaceColor',[0.9 0.9 0.97],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for t = 1:1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    data = avCurves(t).vel;
    plot(plotXAxis,data','Color',color,'LineWidth',lineWidth)
    
    if t > 1
        data = stimLocation(t).optStim(stimLocation(t).optStim > 0)*(70-70/50*t);
        plot(plotXAxis(stimLocation(t).optStim > 0),data,'Color',color,'LineWidth',lineWidth*2)
    end
end
text(0,70,'\fontsize{11}no stim','Color',mymap(1,:))

axis([corridorStart corridorLength 0 70]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.0025,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
rectangle('Position',[-1-viewDistance,0.5,2,80],'FaceColor',[0.9 0.9 0.97],'EdgeColor','none'); % RZ
%rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for t = 1:1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    data = avCurves(t).freq;
    plot(plotXAxis,data,'Color',color,'LineWidth',lineWidth)
    
    if t > 1
        data = stimLocation(t).optStim(stimLocation(t).optStim > 0)*(8-8/50*t);
        plot(plotXAxis(stimLocation(t).optStim > 0),data,'Color',color,'LineWidth',lineWidth*2)
    end
end
text(0,8,'\fontsize{11}no stim','Color',mymap(1,:))


axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick frequency (Hz)');
%legend(plotLegend,'Location','northwest')

%}


% PLOT CURVES
subset = [2, 3, 6];
subset = [2, 3, 4, 5, 6, 7, 8, 9, 10];
subset = [2, 3, 4, 5, 6, 7, 10];
subset = [2, 3, 8, 9];
subset = [2, 3, 4, 6];
subset = [3, 4, 6]; % home box, 60-120, 10-70
subset = [5, 6, 7]; % diff lengths
subset = [4, 6, 10]; % OTP
subset = [4, 5, 6, 7, 10]; % OTP all

rewardZone = sDataFiles{1, 1}.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;
viewDistance = sDataFiles{1, 1}.behavior.viewDistance;  

mymap = [lines(7); 0.5 0.5 0.5; 1 0.1 0.1; 0 1 0.5 ];

plotXAxis = sDataFiles{1, 1}.stats.sessionAvs(1).plotXAxis;


% for plotting stimulus location
for i = 1:1:numel(subset)
    t = subset(i);
    %stimCurves(t) = nanmean(optStimMatrix(sData.trials.trialTypesMeta(nTrialTypes).trials,:));
    stimLocation(t).optStim = avCurves(t).optStim;
    
    stimMax = max(avCurves(t).optStim);
    stimMin = min(avCurves(t).optStim);
    if stimMax < stimMin*1.1 % control
        stimLocation(t).optStim = zeros(size(stimLocation(t).optStim));
    else
        stimLocation(t).optStim(avCurves(t).optStim < (stimMax - stimMin)/3) = 0;
        stimLocation(t).optStim(avCurves(t).optStim >= (stimMax - stimMin)/3) = 1;
    end
end

% stimStart = plotXAxis(min(find(stimLocation == 1)));
% stimEnd = Xax(max(find(stimLocation == 1)));


figure('Color','white','Position',[0 0 850 350]); %pos of figure [left bottom width height]

lineWidth = 1.5;
%suptitle([sData.sessionInfo.sessionID newline])

subplot(1,2,1); % Velocity
hold on

rectangle('Position',[corridorStart+1,0.5-20,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5-20,10,80],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
rectangle('Position',[-1-viewDistance,0.5-20,2,80],'FaceColor',[0.9 0.9 0.97],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
% rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i = 1:1:numel(subset)
    t = subset(i);
    %vr.trialType;
    color = mymap(t,:);
    data = avCurves(t).vel;
    plot(plotXAxis,data','Color',color,'LineWidth',lineWidth)
    
    if t > 1
        data = stimLocation(t).optStim(stimLocation(t).optStim > 0)*(10-30/50*t);
        plot(plotXAxis(stimLocation(t).optStim > 0),data,'Color',color,'LineWidth',lineWidth*2)
    end
end
text(0,70,'\fontsize{11}no stim','Color',mymap(1,:))

axis([corridorStart corridorLength -20 10]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity change (cm/s)');



subplot(1,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.05-3,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.05-3,10,80],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
rectangle('Position',[-1-viewDistance,0.05-3,2,80],'FaceColor',[0.9 0.9 0.97],'EdgeColor','none'); % RZ
%rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i = 1:1:numel(subset)
    t = subset(i);
    
    color = mymap(t,:);
    data = avCurves(t).freq;
    plot(plotXAxis,data,'Color',color,'LineWidth',lineWidth)
    
    if t > 1
        data = stimLocation(t).optStim(stimLocation(t).optStim > 0)*(3-6/50*t);
        plot(plotXAxis(stimLocation(t).optStim > 0),data,'Color',color,'LineWidth',lineWidth*2)
    end
end
text(0,8,'\fontsize{11}no stim','Color',mymap(1,:))


axis([corridorStart corridorLength -3 3]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick rate change (Hz)');
%legend(plotLegend,'Location','northwest')




%% SAVE FIG
% fileName = 'ctrlChangeStandard'; 
% fileName = 'ctrlDiffLengths'; 
% fileName = 'ctrlChangeSummaryOTP1'; 
% fileName = 'ctrlChangeSummaryOTPAll';
% fileName = 'ctrlChangeSummaryHB-OTP1-OTP2'; 
% fileName = 'ctrlChangeSummary1'; 
% fileName = 'ctrlChangeSummary2'; 
% fileName = 'ctrlChangeSummaryAll'; 

% fileName = 'opsinChangeStandard'; 
% fileName = 'opsinChangeDiffLengths'; 

% fileName = 'opsinChangeSummaryOTP1'; 
% fileName = 'opsinChangeSummaryOTPAll';
% fileName = 'opsinChangeSummaryHB-OTP1-OTP2'; 
% fileName = 'opsinChangeSummary1'; 
% fileName = 'opsinChangeSummary2'; 
% fileName = 'opsinChangeSummaryAll'; 
% fileName = 'maskStimChangeSummary'; 


if ~exist([filePath 'summFigs'], 'dir')
    mkdir([filePath 'summFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'summFigs',fileName),'.png'));
close(gcf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate FIRST LICKS (DISTRIBUTION OF FIRST LICKS RELATIVE TO CONTROL)
if exist('curves')
    clear('curves')
end
curves(15).vel = 0;
curves(15).freq = 0;
curves(15).optStim = 0;
curves(15).firstLick = 0;
curves(15).prot = '';

for f = 1:1:nFiles
    
    for i = 1:1:size(sDataFiles{1,f}.trials.trialTypesMeta,2)
        t = sDataFiles{1,f}.trials.trialTypesMeta(i).trialTypeIndicator + 1;
        r = size(curves(t).firstLick,1) + 1;
        
        curves(t).firstLick(r,:) = sDataFiles{1,f}.trials.trialTypesMeta(i).firstLicksCm(2) - sDataFiles{1,f}.trials.trialTypesMeta(1).firstLicksCm(2);
        %curves(t).freq(r,:) = sDataFiles{1,f}.stats.sessionMedians(i).medLickFreqInBin - sDataFiles{1,f}.stats.sessionMedians(1).medLickFreqInBin;
        curves(t).optStim(r,:) = sDataFiles{1,f}.stats.sessionAvs(i).avOptStimMatrix;
        curves(t).prot{r,:} = sDataFiles{1,f}.trials.trialTypesMeta(i).stimProtocol;
        
    end
    
end


% ksdensity(data)

%% EXAMPLE
% data = curves(2).firstLick;
% x_values = -50:5:50;
% pd = fitdist(data,'Kernel'); % or e.g. 'Normal'
% y = pdf(pd,x_values');
% plot(x_values,y,'LineWidth',2)
% equivalent with:
% ksdensity(data)

%{
figure
hold on
x = -150:1:100;
data = curves(2).firstLick;

pd = fitdist(data,'Kernel');
y = pdf(fitdist(data,'Kernel'),x');
plot(x,y,'LineWidth',2)

plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',2)
%}

mymap = [lines(7); 0.5 0.5 0.5; 1 0.1 0.1; 0 1 0.5 ];

figure('Color','white','Position',[0 0 425 350]);
hold on
x = -150:1:100;
lineWidth = 2;

rectangle('Position',[-0.5,0.0005,1,0.1],'FaceColor',[0.1 0.1 0.1],'EdgeColor','none'); % home-box


% data = curves(1).firstLick;
% plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(1,:))

data = curves(2).firstLick;
plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(2,:))
data = curves(3).firstLick;
plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(3,:))
% data = curves(4).firstLick;
% plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(4,:))
% data = curves(5).firstLick;
% plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(5,:))
data = curves(6).firstLick;
plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(6,:))
data = curves(7).firstLick;
plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(7,:))
data = curves(8).firstLick;
plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(8,:))
data = curves(9).firstLick;
plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(9,:))
data = curves(10).firstLick;
plot(x,pdf(fitdist(data,'Kernel'),x'),'LineWidth',lineWidth,'color',mymap(10,:))

ylim([0 0.09])
xlim([-70 50])

xlabel('Shift of first lick position (cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Probability density');




% fileName = 'opsinFirstLickShiftDistrAll'; 
% fileName = 'opsinFirstLickShiftDistrStandardAndDiffLengths'; 
% fileName = 'opsinDiffLengths'; 
% fileName = 'opsinFirstLickShiftDistrStandard4mice'; 
% fileName = 'maskStimFirstLickShiftDistr';

% fileName = 'controlFirstLickShiftDistrAll'; 
% fileName = 'controlFirstLickShiftDistrStandardAndDiffLengths'; 

if ~exist([filePath 'summFigs'], 'dir')
    mkdir([filePath 'summFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'summFigs',fileName),'.png'));
close(gcf)










ksdensity(data)

figure; 
hold on
ksdensity(curves(2).firstLick,'LineWidth',1) % Full tr
%plot(fitdist(curves(2).firstLick,'Normal'))
ksdensity(curves(3).firstLick) % HB
ksdensity(curves(4).firstLick) % OP 60-120
ksdensity(curves(5).firstLick) % 10-20
ksdensity(curves(6).firstLick) % 10-70
ksdensity(curves(7).firstLick) % 10-130
ksdensity(curves(8).firstLick) % -80--40
ksdensity(curves(9).firstLick) % -40-10
ksdensity(curves(10).firstLick) % 10-50


figure; 
hold on
histogram(curves(2).firstLick)
ksdensity(curves(2).firstLick)
histogram(curves(3).firstLick)
ksdensity(curves(3).firstLick)

for t = 1:1:10
    if size(curves(t).optStim,1)>0
        avCurves(t).vel = mean(curves(t).vel);
        avCurves(t).freq = mean(curves(t).freq);
        avCurves(t).optStim = mean(curves(t).optStim);
        avCurves(t).prot = curves(t).prot{1,1};
    end
end


