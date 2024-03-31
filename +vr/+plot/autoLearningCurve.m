%% GET DATA FILES
% files = dir(fullfile(selectedPath,'*.mat'));
%
function learningFig = autoLearningCurve(files)


nSessions = size(files,1);

sDataFiles = cell(1,nSessions);

for f = 1:1:size(files,1)
     
   load(fullfile(files(f).folder,files(f).name));
   % insert analysis function here:
   sDataFiles{f} = sData;
   clear('sData');
      
   
end

filePath = files(1).folder;

%% PLOT DATA


for i = 1:1:numel(sDataFiles)
    
    % various trial lengths
    binNumber = sDataFiles{1,i}.behavior.trialMatrices.meta.binNumber;
    binSize = sDataFiles{1,i}.behavior.trialMatrices.meta.binSize;
    j = find(sDataFiles{1,i}.stats.sessionAvs(1).plotXAxis == binSize);    %60/2 + 1; % first bin of the corridor
    k = j + sDataFiles{1,i}.behavior.rewardZone/binSize - 1; % last bin before RZ   
    
    
    
    quantilesSessions(i,:) = sDataFiles{1,i}.trials.trialTypesMeta(1).lickQuartiles;
    if isfield(sDataFiles{1,i}.stats,'hitRates') % old BCN/NBC task
        hitRates(i,:) = sDataFiles{1,i}.stats.hitRates.hitRateNBC;
    elseif isfield(sDataFiles{1,i}.trials,'trialTypesMeta') % new cue task
        hitRates(i,:) = sDataFiles{1,i}.trials.trialTypesMeta(1).hitRate;
    end
   
    
    lickCounts(i,1) = mean(sum(sDataFiles{1,i}.behavior.trialMatrices.licksInBin(:,j:k)'));
    lickCountsSD(i,1) = std(sum(sDataFiles{1,i}.behavior.trialMatrices.licksInBin(:,j:k)'));
    
    nAllTrials(i,1) = size(sDataFiles{1,i}.behavior.trialMatrices.binVel,1);

    peakFreq(i,1) = nanmean(sDataFiles{1,i}.behavior.trialMatrices.lickFreqInBin(:,k));
    peakFreqSD(i,1) = nanstd(sDataFiles{1,i}.behavior.trialMatrices.lickFreqInBin(:,k));

    startLickFreq(i,1) = nanmean(sDataFiles{1,i}.behavior.trialMatrices.lickFreqInBin(:,j));
    startLickFreqSD(i,1) = nanstd(sDataFiles{1,i}.behavior.trialMatrices.lickFreqInBin(:,j));
     
    velChange(i,1) = mean((sDataFiles{1,i}.behavior.trialMatrices.binVel(:,k) - sDataFiles{1,i}.behavior.trialMatrices.binVel(:,j)) ./ sDataFiles{1,i}.behavior.trialMatrices.binVel(:,j))*100;
    velChangeSD(i,1) = std((sDataFiles{1,i}.behavior.trialMatrices.binVel(:,k) - sDataFiles{1,i}.behavior.trialMatrices.binVel(:,j)) ./ sDataFiles{1,i}.behavior.trialMatrices.binVel(:,j))*100;   
    
    splt = strsplit(sDataFiles{1,i}.sessionInfo.date,'.');
    date = [splt{2} '/' splt{3}];
    sessionDates{i} = date;
    %sessionDates{i} = sDataFiles{1,i}.sessionInfo.date;
    
    splt = strsplit(sDataFiles{1,i}.sessionInfo.sessionID,'-'); 
    splt = strsplit(splt{5},'cm');
    corrLengths(i) = str2num(splt{1});
end





learningFig = figure('Color','white','Position',[0 0 1200 1200]);

subplot(4,1,1);
hold on

rectangle('Position',[-1,79.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

x = 1:1:nSessions;
y25 = quantilesSessions(:,1);
y50 = quantilesSessions(:,2);
y75 = quantilesSessions(:,3);
yneg = y50 - y25;
ypos = y75 - y50;

xticklabels =[{''} sessionDates {''}]; % sessionDates; 
xticks = 0:1:(nSessions + 1); %linspace((1, size(BinVel, 2), numel(xticklabels));

%Plot:
errorbar(x,y50,yneg,ypos,'o-');

axis([0 (nSessions + 1) 0 100]);


set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel('Lick distribution (% of path)');
set(gca, 'TickDir', 'out');
title(['\fontsize{13}' sDataFiles{1,1}.mouseInfo.name ' - Lick distribution']);
%xticklabels = {'light off', 'light on'};
%xticks = 1:1:2;
%set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);



subplot(4,1,3); % HIT RATES
hold on

rectangle('Position',[-1,79.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
plot(hitRates,'o-')

axis([0 (nSessions + 1) 0 110]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);
ylabel('Hit trials (%, reward in reward zone)');
title('\fontsize{13}Hit rate');



subplot(4,1,4); % TRIALS COMPLETED

hold on
rectangle('Position',[-1,79.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

    x = 1:nSessions;
    y = nAllTrials;

plot(x,y,'o-');

    
axis([0 (nSessions + 1) 0 410]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);
title('\fontsize{13}Trials completed');
ylabel('\fontsize{12}Trials completed');
%xlabel('\fontsize{12}Training days');



subplot(4,1,2); % CORRIDOR LENGTHS 
hold on

%rectangle('Position',[-1,79.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

    x = 1:nSessions;
    y = startLickFreq;
    yneg = startLickFreqSD;
    ypos = yneg;    
    errorbar(x,y,yneg,ypos,'o-');

    x = 1:nSessions;
    y = peakFreq;
    yneg = peakFreqSD;
    ypos = yneg;    
    errorbar(x,y,yneg,ypos,'o-');
    
axis([0 (nSessions + 1) 0 8]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);
title('\fontsize{13}Lick freqency at start of the path (blue) and before reward (orange)');
ylabel('\fontsize{11}Lick frequency (Hz)');
%xlabel('\fontsize{12}Training days');
%{
    x = 1:nSessions;
    y = corrLengths;
    
plot(x,y,'o-');

    
axis([0 (nSessions + 1) 0 210]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);
title('\fontsize{13}Distance to reward zone (cm)');
ylabel('\fontsize{12}Distance to reward zone (cm)');
%xlabel('\fontsize{12}Training days');
%}


%{
subplot(3,1,3);
% LICK COUNTS
%plot(lickCounts)

    x = 1:nSessions;
    y = lickCounts;
    yneg = lickCountsSD;
    ypos = yneg;
    
    errorbar(x,y,yneg,ypos,'o-');

    
axis([0 (nSessions + 1) 0 50]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);
title('\fontsize{13}Lick count on the path');
ylabel('\fontsize{12}Lick count (lick/trial)');
%xlabel('\fontsize{12}Training days');
%}

if ~exist([filePath '\summary'], 'dir')
    mkdir([filePath '\summary']);
end

%mkdir([filePath 'summary\']);
saveas(gcf,[filePath '\summary\' sDataFiles{1,1}.mouseInfo.name '-trainingCurves' '.png']);

end
