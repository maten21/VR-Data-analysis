% Plot opto summary for pooled data exps 2020.03.
%
%
%
%
%
%
%
%
%
clear;

[sDataFiles,filePath] = vr.loadData('light');

nFiles = size(sDataFiles,2);


% xticklabels = {'','no stimulus','full trial','Home box','On-the-path',''};
% xticklabels = {'','no stimulus','10 - 20 cm','10 - 70 cm','10 -130 cm',''};
% xticklabels = {'','no stimulus','60 - 120 cm',''};
% xticklabels = {'','no stimulus','10 - 70 cm',''};
% xticklabels = {'','no stimulus','10 - 70 cm','no stimulus','60 - 120 cm',''};
% xticklabels = {'','no stimulus','-80 - 40','-30 - 10 cm','10 - 50 cm',''};
%% Extract data
nTrialTypes = size(sDataFiles{1, 1}.trials.trialTypesMeta,2);
firstTrackBin = 71;
lastTrackBin = 160;
for f = 1:1:nFiles
    
    trFreqs = nanmean(sDataFiles{1, f}.behavior.trialMatrices.lickFreqInBin(:,firstTrackBin:lastTrackBin)');
    trVels = nanmean(sDataFiles{1, f}.behavior.trialMatrices.binVel(:,firstTrackBin:lastTrackBin)');
    trTimes = nansum(sDataFiles{1, f}.behavior.trialMatrices.timeInBin(:,:)');
    for t = 1:1:nTrialTypes
        trials = sDataFiles{1, f}.trials.trialTypesMeta(t).trials;
        
        distr(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).lickQuartiles *1.8;
        firstL(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).firstLicksCm;
        %acc(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).accBeforeRZ;
        hits(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).hitRate;
        
        freq(t,:) = nanmean(trFreqs(trials));
        freqStd(t,:) = std(trFreqs(trials));
        vel(t,:) =  nanmean(trVels(trials));
        velStd(t,:) =  std(trVels(trials));
        timeInBins(t,:) = nanmean(trTimes(trials));
        timeInBinsStd(t,:) = std(trTimes(trials));
    end
    
    allSessions.mouseNames{f} = sDataFiles{1, f}.mouseInfo.name;
    allSessions.lickDistr25(f,:) = distr(:,1);
    allSessions.lickDistr50(f,:) = distr(:,2);
    allSessions.lickDistr75(f,:) = distr(:,3);
    allSessions.firstLicks25(f,:) = firstL(:,1);
    allSessions.firstLicks50(f,:) = firstL(:,2);
    allSessions.firstLicks75(f,:) = firstL(:,3);
    %allSessions.acceleration{f} = acc;
    allSessions.hitRate(f,:) = hits;
    
    allSessions.trialFreqs(f,:) = freq;
    allSessions.trialFreqsStd(f,:) = freqStd;
    allSessions.trialVels(f,:) = vel;
    allSessions.trialVelsStd(f,:) = velStd;
    allSessions.trialTimes(f,:) = timeInBins;
    allSessions.trialTimesStd(f,:) = timeInBinsStd;
    
    allSessions.OTPProt{f,:} = sDataFiles{1, f}.trials.trialTypesMeta(4).stimProtocol;
    
end
clear('distr','firstL','hits','freq','freqStd','vel','velStd','timeInBins','timeInBinsStd')

%% Extract data to compare On-the-path protocols
%{
nTrialTypes = size(sDataFiles{1, 1}.trials.trialTypesMeta,2);
firstTrackBin = 71;
lastTrackBin = 160;
for f = 1:1:nFiles
    
    trFreqs = nanmean(sDataFiles{1, f}.behavior.trialMatrices.lickFreqInBin(:,firstTrackBin:lastTrackBin)');
    trVels = nanmean(sDataFiles{1, f}.behavior.trialMatrices.binVel(:,firstTrackBin:lastTrackBin)');
    trTimes = nansum(sDataFiles{1, f}.behavior.trialMatrices.timeInBin(:,:)');
    for t = 1:1:nTrialTypes
            
        if t == 1
        trials = sDataFiles{1, f}.trials.trialTypesMeta(t).trials;
        
        distr(1,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).lickQuartiles *1.8; % convert to cm
        firstL(1,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).firstLicksCm;
        hits(1,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).hitRate;
        
        freq(1,:) = nanmean(trFreqs(trials));
        freqStd(1,:) = std(trFreqs(trials));
        vel(1,:) =  nanmean(trVels(trials));
        velStd(1,:) =  std(trVels(trials));
        timeInBins(1,:) = nanmean(trTimes(trials));
        timeInBinsStd(1,:) = std(trTimes(trials));
                
        elseif strcmp(sDataFiles{1, f}.trials.trialTypesMeta(t).stimProtocol,'20 Hz sinus from 10 to 70 int. 100 % duty cycle 0 %')
        trials = sDataFiles{1, f}.trials.trialTypesMeta(t).trials;
        
        distr(2,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).lickQuartiles *1.8; % convert to cm
        firstL(2,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).firstLicksCm;
        hits(2,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).hitRate;
        
        freq(2,:) = nanmean(trFreqs(trials));
        freqStd(2,:) = std(trFreqs(trials));
        vel(2,:) =  nanmean(trVels(trials));
        velStd(2,:) =  std(trVels(trials));
        timeInBins(2,:) = nanmean(trTimes(trials));
        timeInBinsStd(2,:) = std(trTimes(trials));
        end
        
    end
    
    allSessions.mouseNames{f} = sDataFiles{1, f}.mouseInfo.name;
    allSessions.lickDistr25(f,:) = distr(:,1);
    allSessions.lickDistr50(f,:) = distr(:,2);
    allSessions.lickDistr75(f,:) = distr(:,3);
    allSessions.firstLicks25(f,:) = firstL(:,1);
    allSessions.firstLicks50(f,:) = firstL(:,2);
    allSessions.firstLicks75(f,:) = firstL(:,3);
    %allSessions.acceleration{f} = acc;
    allSessions.hitRate(f,:) = hits;
    
    allSessions.trialFreqs(f,:) = freq;
    allSessions.trialFreqsStd(f,:) = freqStd;
    allSessions.trialVels(f,:) = vel;
    allSessions.trialVelsStd(f,:) = velStd;
    allSessions.trialTimes(f,:) = timeInBins;
    allSessions.trialTimesStd(f,:) = timeInBinsStd;
    
    allSessions.OTPProt{f,:} = sDataFiles{1, f}.trials.trialTypesMeta(4).stimProtocol;
    
end
clear('distr','firstL','hits','freq','freqStd','vel','velStd','timeInBins','timeInBinsStd')

nTrialTypes = 2;
%}



%% Average data from the same mice

g = 1;
h = 1;
for f = 1:1:nFiles %calculate average matrices
    
    if f == nFiles || ~strcmp(allSessions.mouseNames{f+1},allSessions.mouseNames{g})
    
    distr25 = zeros(size(allSessions.lickDistr50(1,:)));
    distr50 = zeros(size(allSessions.lickDistr50(1,:)));
    distr75 = zeros(size(allSessions.lickDistr50(1,:)));
    firstL25 = zeros(size(allSessions.lickDistr50(1,:)));
    firstL50 = zeros(size(allSessions.lickDistr50(1,:)));
    firstL75 = zeros(size(allSessions.lickDistr50(1,:)));
    hits = zeros(size(allSessions.lickDistr50(1,:)));
    freqs = zeros(size(allSessions.lickDistr50(1,:)));
    freqsStd = zeros(size(allSessions.lickDistr50(1,:)));
    vels = zeros(size(allSessions.lickDistr50(1,:)));
    velsStd = zeros(size(allSessions.lickDistr50(1,:)));
    times = zeros(size(allSessions.lickDistr50(1,:)));
    timesStd = zeros(size(allSessions.lickDistr50(1,:)));
    for i = 1:1:f-g+1
        distr25 = distr25 + allSessions.lickDistr25(g+i-1,:);
        distr50 = distr50 + allSessions.lickDistr50(g+i-1,:);
        distr75 = distr75 + allSessions.lickDistr75(g+i-1,:);
        firstL25 = firstL25 + allSessions.firstLicks25(g+i-1,:);
        firstL50 = firstL50 + allSessions.firstLicks50(g+i-1,:);
        firstL75 = firstL75 + allSessions.firstLicks75(g+i-1,:);
        hits = hits + allSessions.hitRate(g+i-1,:);
        freqs = freqs + allSessions.trialFreqs(g+i-1,:);
        freqsStd = freqsStd + allSessions.trialFreqsStd(g+i-1,:);
        vels = vels + allSessions.trialVels(g+i-1,:);
        velsStd = velsStd + allSessions.trialVelsStd(g+i-1,:);
        times = times + allSessions.trialTimes(g+i-1,:);
        timesStd = timesStd + allSessions.trialTimesStd(g+i-1,:);
    end
       
    
    mouseAvs.lickDistr25(h,:) = distr25/i;
    mouseAvs.lickDistr50(h,:) = distr50/i;
    mouseAvs.lickDistr75(h,:) = distr75/i;
    mouseAvs.firstLicks25(h,:) = firstL25/i;
    mouseAvs.firstLicks50(h,:) = firstL50/i;
    mouseAvs.firstLicks75(h,:) = firstL75/i;
    
    mouseAvs.hitRate(h,:) = hits/i;
    mouseAvs.trialFreqs(h,:) = freqs/i;
    mouseAvs.trialFreqsStd(h,:) = freqsStd/i;
    mouseAvs.trialVels(h,:) = vels/i;
    mouseAvs.trialVelsStd(h,:) = velsStd/i;
    mouseAvs.trialTimes(h,:) = times/i;
    mouseAvs.trialTimesStd(h,:) = timesStd/i;
    
    h = h+1;
    g = f+1;
    end
end
clear('distr25','distr50','distr75','firstL25','firstL50','firstL75','hits','freqs','freqsStd','vels','velsStd','times','timesStd')

trialType = 2;
for i =1:1:size(mouseAvs.firstLicks50,1)
    mouseAvsNormEffect.lickDistr25(i,1) = mouseAvs.lickDistr25(i,trialType)/mouseAvs.lickDistr25(i,1);
    mouseAvsNormEffect.lickDistr50(i,1) = mouseAvs.lickDistr50(i,trialType)/mouseAvs.lickDistr50(i,1);
    mouseAvsNormEffect.lickDistr75(i,1) = mouseAvs.lickDistr75(i,trialType)/mouseAvs.lickDistr75(i,1);
    mouseAvsNormEffect.firstLicks25(i,1) = mouseAvs.firstLicks25(i,trialType)/mouseAvs.firstLicks25(i,1);
    mouseAvsNormEffect.firstLicks50(i,1) = mouseAvs.firstLicks50(i,trialType)/mouseAvs.firstLicks50(i,1);
    mouseAvsNormEffect.firstLicks75(i,1) = mouseAvs.firstLicks75(i,trialType)/mouseAvs.firstLicks75(i,1);
    
    mouseAvsNormEffect.hitRate(i,1) = mouseAvs.hitRate(i,trialType)/mouseAvs.hitRate(i,1);
    mouseAvsNormEffect.trialFreqs(i,1) = mouseAvs.trialFreqs(i,trialType)/mouseAvs.trialFreqs(i,1);
    mouseAvsNormEffect.trialFreqsStd(i,1) = mouseAvs.trialFreqsStd(i,trialType)/mouseAvs.trialFreqsStd(i,1);
    mouseAvsNormEffect.trialVels(i,1) = mouseAvs.trialVels(i,trialType)/mouseAvs.trialVels(i,1);
    mouseAvsNormEffect.trialVelsStd(i,1) = mouseAvs.trialVelsStd(i,trialType)/mouseAvs.trialVelsStd(i,1);
    mouseAvsNormEffect.trialTimes(i,1) = mouseAvs.trialTimes(i,trialType)/mouseAvs.trialTimes(i,1);
    mouseAvsNormEffect.trialTimesStd(i,1) = mouseAvs.trialTimesStd(i,trialType)/mouseAvs.trialTimesStd(i,1);
end




%%
% allSessions = allSessions60; Xax = [3 4];
% allSessions = allSessions10; Xax = [1 2];
% Xax = [1 2];


%% Plot fig INDIVIDUAL SESSIONS
Xax = 1:1:nTrialTypes;
color = lines(1);
% Xcount = nTrialTypes; Xcount = 4;
xticks = 0:1:(Xcount + 1);

figure('Color','white','Position',[0 0 900 600]);

%% Lick distribution
subplot(2,3,1) % Lick distribution

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y25 =  allSessions.lickDistr25(f,:);
    y50 =  allSessions.lickDistr50(f,:);
    y75 =  allSessions.lickDistr75(f,:);
    yneg = y50 - y25;
    ypos = y75 - y50;
    
errorbar(Xax,y50,yneg,ypos,'o-','color',color);

end

axis([0 (Xcount + 1) 0 180]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Lick distr. on track(cm, quartiles)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% First licks
subplot(2,3,2) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y25 = allSessions.firstLicks25(f,:);
    y50 = allSessions.firstLicks50(f,:);
    y75 = allSessions.firstLicks75(f,:);
    yneg = y50 - y25;
    ypos = y75 - y50;
    
errorbar(Xax,y50,yneg,ypos,'o-','color',color);

end

axis([0 (Xcount + 1) 0 180]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Fitst lick on track (cm, quartiles)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Hit Rate
subplot(2,3,3) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles
   
    plot(Xax,allSessions.hitRate(f,:),'o-','color',color);

end

axis([0 (Xcount + 1) 40 105]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Hit rate (%)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Average velocity
subplot(2,3,4) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y = allSessions.trialVels(f,:);
    yStd = allSessions.trialVelsStd(f,:);

    errorbar(Xax,y,yStd,'o-','color',color);

end

axis([0 (Xcount + 1) 0 80]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Velocity on track (cm/s, mean ± SD )'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Average lick frequency
subplot(2,3,5) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y = allSessions.trialFreqs(f,:);
    yStd = allSessions.trialFreqsStd(f,:);

    errorbar(Xax,y,yStd,'o-','color',color);

end

axis([0 (Xcount + 1) 0 5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Lick rate on track (Hz, mean ± SD )'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Trial time 
subplot(2,3,6) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y = allSessions.trialTimes(f,:);
    yStd = allSessions.trialTimesStd(f,:);

    
    errorbar(Xax,y,yStd,'o-','color',color);

end

axis([0 (Xcount + 1) 0 30]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Trial time (s, mean ± SD )'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');



%%
% mouseAvs = mouseAvs60; Xax = [3 4];
% mouseAvs = mouseAvs10; Xax = [1 2];
% xticklabels = {'','10 - 70 cm','60 - 120 cm',''};
% Xcount = 2;
% xticks = 0:1:(Xcount + 1);
% Xax = [1 2];

% mouseAvs = mouseAvsNormEffect10;
% mouseAvs = mouseAvsNormEffect60;
% yData1 = mouseAvsNormEffect10;
% yData2 = mouseAvsNormEffect60;
%% Plot fig AVERAGE SESSIONS

figure('Color','white','Position',[0 0 900 600]);

Xax = 1:1:nTrialTypes;
color = lines(1);
% Xcount = nTrialTypes; Xcount = 4;
xticks = 0:1:(Xcount + 1);
nMice = size(mouseAvs.firstLicks50,1);
%% Lick distribution
subplot(2,3,1) % Lick distribution

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i =1:1:nMice

    y25 = mouseAvs.lickDistr25(i,:);
    y50 = mouseAvs.lickDistr50(i,:);
    y75 = mouseAvs.lickDistr75(i,:);
    yneg = y50 - y25;
    ypos = y75 - y50;
    
    % errorbar(Xax,y50,yneg,ypos,'o-','color',color);
    plot(Xax,y50,'o-','color',color);
    
end

% plot(Xax,mean(mouseAvs.lickDistr50),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)
errorbar(Xax,mean(mouseAvs.lickDistr50),std(mouseAvs.lickDistr50),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10,'CapSize',12,'Linewidth',1.5)


axis([0 (Xcount + 1) 0 180]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Lick distr. on track (%, quartiles)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% First licks
subplot(2,3,2) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i =1:1:nMice

    y25 = mouseAvs.firstLicks25(i,:);
    y50 = mouseAvs.firstLicks50(i,:);
    y75 = mouseAvs.firstLicks75(i,:);
    yneg = y50 - y25;
    ypos = y75 - y50;
    
    % errorbar(Xax,y50,yneg,ypos,'o-','color',color);
    plot(Xax,y50,'o-','color',color);
end

% plot(Xax,mean(mouseAvs.firstLicks50),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)
errorbar(Xax,mean(mouseAvs.firstLicks50),std(mouseAvs.firstLicks50),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10,'CapSize',12,'Linewidth',1)

axis([0 (Xcount + 1) 0 180]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Fitst lick on track (cm, quartiles)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Hit Rate
subplot(2,3,3) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i =1:1:nMice
   
    plot(Xax,mouseAvs.hitRate(i,:),'o-','color',color);

end

% plot(Xax,mean(mouseAvs.hitRate),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)
errorbar(Xax,mean(mouseAvs.hitRate),std(mouseAvs.hitRate),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10,'CapSize',12,'Linewidth',1)

axis([0 (Xcount + 1) 40 105]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Hit rate (%)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Average velocity

subplot(2,3,4) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i =1:1:nMice

    y = mouseAvs.trialVels(i,:);
    yStd = mouseAvs.trialVelsStd(i,:);
    
    % errorbar(Xax,y,yStd,'o-','color',color);
    plot(Xax,y,'o-','color',color);
end
%plot(Xax,mean(mouseAvs.trialVels),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)
errorbar(Xax,mean(mouseAvs.trialVels),std(mouseAvs.trialVels),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10,'CapSize',12,'Linewidth',1)

axis([0 (Xcount + 1) 0 80]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Velocity on track (cm/s, mean ± SD )'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Average lick frequency

subplot(2,3,5) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i =1:1:nMice

    y = mouseAvs.trialFreqs(i,:);
    yStd = mouseAvs.trialFreqsStd(i,:);
    
    % errorbar(Xax,y,yStd,'o-','color',color);
    plot(Xax,y,'o-','color',color);
end
%plot(Xax,mean(mouseAvs.trialFreqs),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)
errorbar(Xax,mean(mouseAvs.trialFreqs),std(mouseAvs.trialFreqs),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10,'CapSize',12,'Linewidth',1)

axis([0 (Xcount + 1) 0 5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Lick rate on track (Hz, mean ± SD )'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Trial time 

subplot(2,3,6) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i =1:1:nMice
    
    y = mouseAvs.trialTimes(i,:);
    yStd = mouseAvs.trialTimesStd(i,:);
    
    % errorbar(Xax,y,yStd,'o-','color',color);
    plot(Xax,y,'o-','color',color);
end
%plot(Xax,mean(mouseAvs.trialTimes),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)
errorbar(Xax,mean(mouseAvs.trialTimes),std(mouseAvs.trialTimes),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10,'CapSize',12,'Linewidth',1)

axis([0 (Xcount + 1) 0 30]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Trial time (s, mean ± SD )'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');


saveas(gcf,strcat(fullfile(filePath,'summFigs','mouseAverageComp'),'.png'));



%% Correlation

%{
figure
hold on
for f =1:1:nFiles
   
    plot((firstLicks{f}(4,2)-firstLicks{f}(1,2)),hitRate{f}(4),'o','color',color);

end



trType = 1;
xMin = min(data1(:,trType))*0.95;
xMax = max(data1(:,trType))*1.05;
yMin = 140;
yMax = 180;


figure('Color','white','Position',[0 0 300 300]); 
hold on
plot(data1(:,trType),data2(:,trType),'o')

[fitX, fitY, R2, P] = linFit(data1(:,trType),data2(:,trType));
if P < 0.05
    plot(fitX,fitY,'-r')
    text(xMin+(xMax-xMin)/3,yMin+(yMax-yMin)/3,['R^2 = ' num2str(R2) newline 'P = ' num2str(P)])
end
ylim([yMin yMax])
xlim([xMin xMax])
title(xticklabels{trType+1})

% xlabel('Lick distr. on track (cm, median)')
% ylabel('Hit rate (%)')
% ylabel('Lick rate (Hz)')
% ylabel('Lick distr. on track (cm)')
% ylabel('First lick on track (cm)')
 xlabel('Velocity (cm/s)')
% xlabel('Mean trial duration(s)')
% xlabel('Mean trial duration(s)')
% xlabel('Lick rate (Hz)')
% xlabel('First lick on track (cm)')

%}

% data1 = allSessions.trialTimes./allSessions.hitRate;


figure('Color','white','Position',[0 0 900 300]);


mymap = lines(1);
color = mymap(1,:); 

% data1 = allSessions.firstLicks50; Xlabel = 'First lick in control (cm, median)';
% data1 = allSessions.lickDistr50; Xlabel = 'Lick distribution in control (cm, median)';
% data1 = allSessions.hitRate; Xlabel = 'Hitrate in control (%)';
% data1 = allSessions.trialVels; Xlabel = 'Velocity in control (cm/s)';
% data1 = allSessions.trialFreqs; Xlabel = 'Lick rate in control (Hz)';
% data1 = allSessions.trialTimes; Xlabel = 'Trial time in control (s)';

% data1 = allSessions.trialVels; Xlabel = 'Change in velocity (cm/s)';
% data1 = allSessions.trialFreqs; Xlabel = 'Change in lick rate (Hz)';

% data2 = allSessions.firstLicks50; Ylabel = 'Shift of first lick';
% data2 = allSessions.lickDistr50; Ylabel = 'Shift of lick distr.';
% data2 = allSessions.hitRate; Ylabel = 'Change in hit rate';
% data2 = allSessions.trialVels; Ylabel = 'Change in velocity (cm/s)';
% data2 = allSessions.trialFreqs; Ylabel = 'Change in lick rate (Hz)';
% data2 = allSessions.trialTimes; Ylabel = 'Change in trial duration (s)';

 xMin = -1;
 xMax = 2;
% xMin = min(data1(:,1))*0.95;
% xMax = max(data1(:,1))*1.05;
yMin = -30;
yMax = 35;
subplot(1,3,1)
hold on

%x = data1(:,1); 
x = (data1(:,2)-data1(:,1));
y = (data2(:,2)-data2(:,1));
plot(x,y,'o','color',color);
[fitX, fitY, R2, P] = linFit(x,y);
if P < 0.05
    plot(fitX,fitY,'-r')
    text(xMin+(xMax-xMin)/3,yMin+(yMax-yMin)/3,['R^2 = ' num2str(R2) newline 'P = ' num2str(P)])
end

ylim([yMin yMax])
xlim([xMin xMax])
title(xticklabels{3})
ylabel(Ylabel)

subplot(1,3,2)
hold on

% x = data1(:,1);
x = (data1(:,3)-data1(:,1));
y = (data2(:,3)-data2(:,1));
plot(x,y,'o','color',color);
[fitX, fitY, R2, P] = linFit(x,y);
if P < 0.05
    plot(fitX,fitY,'-r')
    text(xMin+(xMax-xMin)/2,yMin+(yMax-yMin)/2,['R^2 = ' num2str(R2) newline 'P = ' num2str(P)])
end

ylim([yMin yMax])
xlim([xMin xMax])
title(xticklabels{4})
xlabel(Xlabel)

subplot(1,3,3)
hold on

%x = data1(:,1);
x = (data1(:,4)-data1(:,1));
y = (data2(:,4)-data2(:,1));
plot(x,y,'o','color',color);
[fitX, fitY, R2, P] = linFit(x,y);
if P < 0.05
    plot(fitX,fitY,'-r')
    text(xMin+(xMax-xMin)/2,yMin+(yMax-yMin)/2,['R^2 = ' num2str(R2) newline 'P = ' num2str(P)])
end

ylim([yMin yMax])
xlim([xMin xMax])
title(xticklabels{5})


%
%%%
%%%%%
%%%

%% Extract data for statistics save excell
excellFileName = ['C:\Users\Mate Neubrandt\Dropbox (UIO Physiology Dropbox)\MateData\RECORDINGS\exp Opto 2020 winter\' 'summaryTable'];
% sheet = 'ctrlStandard';
% sheet = 'ctrlDiffLengths';
% sheet = 'opsinStandard';
% sheet = 'opsinDiffLengths';

header = xticklabels(2:size(xticklabels,2)-1);

name = {'Lick distr.'}; 
data = mouseAvs.lickDistr50;
counter = 0;
range = ['A' num2str(counter+1)]; 
xlswrite(excellFileName,name,sheet,range)
range = ['A' num2str(counter+2) ':D' num2str(counter+2)];
xlswrite(excellFileName,header,sheet,range)
range = ['A' num2str(counter+3) ':D' num2str(counter+3 + size(data,1)-1)];
xlswrite(excellFileName,data,sheet,range)

name = {'First lick'}; 
data = mouseAvs.firstLicks50;
counter = counter + 3 + size(data,1);
range = ['A' num2str(counter+1)]; 
xlswrite(excellFileName,name,sheet,range)
range = ['A' num2str(counter+2) ':D' num2str(counter+2)];
xlswrite(excellFileName,header,sheet,range)
range = ['A' num2str(counter+3) ':D' num2str(counter+3 + size(data,1)-1)];
xlswrite(excellFileName,data,sheet,range)

name = {'Hit rate'}; 
data = mouseAvs.hitRate;
counter = counter + 3 + size(data,1);
range = ['A' num2str(counter+1)]; 
xlswrite(excellFileName,name,sheet,range)
range = ['A' num2str(counter+2) ':D' num2str(counter+2)];
xlswrite(excellFileName,header,sheet,range)
range = ['A' num2str(counter+3) ':D' num2str(counter+3 + size(data,1)-1)];
xlswrite(excellFileName,data,sheet,range)

name = {'Velocity'}; 
data = mouseAvs.trialVels;
counter = counter + 3 + size(data,1);
range = ['A' num2str(counter+1)]; 
xlswrite(excellFileName,name,sheet,range)
range = ['A' num2str(counter+2) ':D' num2str(counter+2)];
xlswrite(excellFileName,header,sheet,range)
range = ['A' num2str(counter+3) ':D' num2str(counter+3 + size(data,1)-1)];
xlswrite(excellFileName,data,sheet,range)

name = {'Lick rate'}; 
data = mouseAvs.trialFreqs;
counter = counter + 3 + size(data,1);
range = ['A' num2str(counter+1)]; 
xlswrite(excellFileName,name,sheet,range)
range = ['A' num2str(counter+2) ':D' num2str(counter+2)];
xlswrite(excellFileName,header,sheet,range)
range = ['A' num2str(counter+3) ':D' num2str(counter+3 + size(data,1)-1)];
xlswrite(excellFileName,data,sheet,range)

name = {'Trial times'}; 
data = mouseAvs.trialTimes;
counter = counter + 3 + size(data,1);
range = ['A' num2str(counter+1)]; 
xlswrite(excellFileName,name,sheet,range)
range = ['A' num2str(counter+2) ':D' num2str(counter+2)];
xlswrite(excellFileName,header,sheet,range)
range = ['A' num2str(counter+3) ':D' num2str(counter+3 + size(data,1)-1)];
xlswrite(excellFileName,data,sheet,range)





%%%%%%
%%%%%%%%%%%
%%%%%%%%%%%%%%%%% SEE:  vr.plot.OptoSessionCompSummaryCurves !!! that is
%%%%%%%%%%%%%%%%% the latest!
%%%%%%%%%%%
%%%%%%

%% Generate summary curves ALL sessions





for t = 1:1:nTrialTypes
    
    tempVel = 0;
    tempLickFreq = 0;
    tempOptStim = 0;
    for f = 1:1:nFiles
        
        tempVel = tempVel + sDataFiles{1, f}.stats.sessionMedians(t).medBinVel;
        tempLickFreq = tempLickFreq + sDataFiles{1, f}.stats.sessionMedians(t).medLickFreqInBin;
        tempOptStim = tempOptStim + sDataFiles{1, f}.stats.sessionAvs(t).avOptStimMatrix; 
    end
    avCurves(t).vel = tempVel / nFiles;
    avCurves(t).freq = tempLickFreq / nFiles;
    avCurves(t).optStim = tempOptStim / nFiles;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correct trial indicator errors!!!
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
subset = [1, 2, 3, 4, 5, 6, 7, 10];
subset = [1, 2, 3, 8, 9];
subset = [1, 2, 3, 4, 6];
subset = [1, 5, 6, 7]; % diff lengths

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

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
rectangle('Position',[-1-viewDistance,0.5,2,80],'FaceColor',[0.9 0.9 0.97],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i = 1:1:numel(subset)
    t = subset(i);
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

rectangle('Position',[corridorStart+1,0.05,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.05,10,80],'FaceColor',[0.9 0.97 0.9],'EdgeColor','none'); % RZ
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
% fileName = 'ctrlStandard2'; 
% fileName = 'ctrlDiffLengths'; 
% fileName = 'controlSummary1'; 
% fileName = 'controlSummary2'; 
% fileName = 'opsinStandard'; 
% fileName = 'opsinDiffLengths'; 
% fileName = 'opsinShortStim'; 
% fileName = 'opsinSummary1'; 
% fileName = 'opsinSummary2'; 


if ~exist([filePath 'summFigs'], 'dir')
    mkdir([filePath 'summFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'summFigs',fileName),'.png'));
close(gcf)








figure
hold on

plot(avCurves(1).vel)
plot(avCurvesOpsinStandard(1).vel)
plot(avCurvesOpsinDiffLengths(1).vel)

plot(avCurvesOpsinStandard(3).vel)
plot(avCurvesOpsinStandard(4).vel)


figure
hold on
plot(avCurves(1).freq)
plot(avCurvesOpsinStandard(1).freq)
plot(avCurvesOpsinDiffLengths(1).freq)

plot(avCurvesOpsinStandard(3).freq)
plot(avCurvesOpsinStandard(4).freq)























%% Reward time

rewardTimes =  allSessions.trialTimes./(allSessions.hitRate/100);


%% Plot
figure('Color','white','Position',[0 0 300 300]);

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for i =1:1:nMice
   
    plot(Xax,rewardTimes(i,:),'o-','color',color);

end

plot(Xax,mean(rewardTimes),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)

axis([0 (nTrialTypes + 1) 0 40]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Task efficacy (s/reward)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');








%% Plot EFFECT AVERAGE SESSIONS %

figure('Color','white','Position',[0 0 900 600]);

nMice = size(mouseAvs.firstLicks50,1);
%% Lick distribution
subplot(2,3,1) % Lick distribution

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
y = [yData1.lickDistr50 yData2.lickDistr50];
for i =1:1:nMice
   
    plot(Xax,y(i,:),'o-','color',color);

end

plot(Xax,mean(y),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)

axis([0 (Xcount + 1) 0.5 1.5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Lick distr. on track (relative to control)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% First licks
subplot(2,3,2) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
y = [yData1.firstLicks50(:,:) yData2.firstLicks50(:,:)];
for i =1:1:nMice
    
    plot(Xax,y(i,:),'o-','color',color);

end

plot(Xax,mean(y),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)

axis([0 (Xcount + 1) 0.5 1.5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Fitst lick on track (relative to control)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Hit Rate
subplot(2,3,3) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
y = [yData1.hitRate(:,:) yData2.hitRate(:,:)];
for i =1:1:nMice
   
    plot(Xax,y(i,:),'o-','color',color);

end

plot(Xax,mean(y),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)

axis([0 (Xcount + 1) 0.5 1.5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Hit rate (relative to control)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Average velocity

subplot(2,3,4) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
y = [yData1.trialVels(:,:) yData2.trialVels(:,:)];
for i =1:1:nMice
    
    plot(Xax,y(i,:),'o-','color',color);

end
plot(Xax,mean(y),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)

axis([0 (Xcount + 1) 0.5 1.5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Velocity on track (relative to control)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Average lick frequency

subplot(2,3,5) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
y = [yData1.trialFreqs(:,:) yData2.trialFreqs(:,:)];
for i =1:1:nMice
    
    plot(Xax,y(i,:),'o-','color',color);

end
plot(Xax,mean(y),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)

axis([0 (Xcount + 1) 0.5 1.5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Lick rate on track (relative to control)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Trial time 

subplot(2,3,6) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
y = [yData1.trialTimes(:,:) yData2.trialTimes(:,:)];
for i =1:1:nMice
    
    plot(Xax,y(i,:),'o-','color',color);

end
plot(Xax,mean(y),'o-','color',lines(1),'markerfacecolor',lines(1),'MarkerSize',10)

axis([0 (Xcount + 1) 0.5 1.5]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Trial time (relative to control)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');




