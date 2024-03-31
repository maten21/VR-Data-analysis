% Plot opto summary for pooled data exps 2019.11.
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

sDataFiles = vr.loadData('light');

nFiles = size(sDataFiles,2);


%% Extract data
for f = 1:1:nFiles

    nTrialTypes = size(sDataFiles{1, 1}.trials.trialTypesMeta,2);
    for t = 1:1:nTrialTypes
        distr(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).lickQuartiles;
        firstL(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).firstLicksPerc;
        acc(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).accBeforeRZ;
        hits(t,:) = sDataFiles{1, f}.trials.trialTypesMeta(t).hitRate;
    end
    lickDistr{f} = distr;
    firstLicks{f} = firstL;
    acceleration{f} = acc;
    hitRate{f} = hits;
    
end


%% Plot fig
Xax = 1:1:nTrialTypes;
color = lines(1);
xticks = 0:1:(nTrialTypes + 1);
xticklabels = {'','no stimulus','full trial','0 - 40 cm','80 - 120 cm','160 - 200 cm',''};


figure('Color','white','Position',[0 0 700 700]);

%% Lick distribution
subplot(2,2,1) % Lick distribution

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y25 = lickDistr{f}(:,1);
    y50 = lickDistr{f}(:,2);
    y75 = lickDistr{f}(:,3);
    yneg = y50 - y25;
    ypos = y75 - y50;
    
errorbar(Xax,y50,yneg,ypos,'o-','color',color);

end

axis([0 (nTrialTypes + 1) 50 100]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Lick distribution (%, quartiles)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');


%% First licks
subplot(2,2,2) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y25 = firstLicks{f}(:,1);
    y50 = firstLicks{f}(:,2);
    y75 = firstLicks{f}(:,3);
    yneg = y50 - y25;
    ypos = y75 - y50;
    
errorbar(Xax,y50,yneg,ypos,'o-','color',color);

end

axis([0 (nTrialTypes + 1) 0 100]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Fitst lick on track (%, quartiles)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');

%% Acceleration
subplot(2,2,3) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles

    y25 = acceleration{f}(:,1);
    y50 = acceleration{f}(:,2);
    y75 = acceleration{f}(:,3);
    yneg = y50 - y25;
    ypos = y75 - y50;
    
errorbar(Xax,y50,yneg,ypos,'o-','color',color);

end

axis([0 (nTrialTypes + 1) -40 20]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Acceleration (cm/s^2, quartiles)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');


%% Hit Rate
subplot(2,2,4) 

hold on
%rectangle('Position',[-1,-0.5,100,1],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for f =1:1:nFiles
   
    plot(Xax,hitRate{f},'o-','color',color);

end

axis([0 (nTrialTypes + 1) 40 105]);

set(gca, 'TickDir', 'out','XTick', xticks,'XTickLabel', xticklabels);
xtickangle(45);

ylabel({'\fontsize{11}Hit rate (%)'});
set(gca, 'TickDir', 'out');
%title('\fontsize{13}Acceleration in the last bin before the reward');














