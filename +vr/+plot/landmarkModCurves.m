%% Plot landmark modulation curves 

function [] = landmarkModCurves(fileName,filePath)

load(fullfile(filePath,fileName));

lickFreqInBin = sData.behavior.trialMatrices.lickFreqInBin;
binVel = sData.behavior.trialMatrices.binVel;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;
Xax = sData.stats.sessionAvs(1).plotXAxis;

rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;

nTrialTypes = size(sData.trials.trialTypesMeta,2);



myMap = lines;
smoothSpan = 5;
lineWidth = 1;

    

figure('Color','white','Position',[0 0 850 350*nTrialTypes])
hold on
for t = 1:1:nTrialTypes

    plotLegend = {['from -140 cm (n = ' num2str(numel(sData.trials.trialTypesMeta(t).randStartPos(1).trials  )) ')'];...
        ['from -120 cm (n = ' num2str(numel(sData.trials.trialTypesMeta(t).randStartPos(2).trials  )) ')'];...
        ['from -100 cm (n = ' num2str(numel(sData.trials.trialTypesMeta(t).randStartPos(3).trials  )) ')'];...
        ['from -80 cm (n = ' num2str(numel(sData.trials.trialTypesMeta(t).randStartPos(4).trials  )) ')']};

    
    % for plotting stimulus location
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(t).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimTrial = abs(stimMax-stimMin)>0.1;
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));
    
    

    sp1 = subplot(nTrialTypes,2,t*2-1);
    hold on
    
    rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
    rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
    if stimTrial
    rectangle('Position',[stimStart,76,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
    end
    
    for i = 1:1:4
        color = myMap(i,:);
        plot(Xax,smoothdata(nanmedian(binVel(sData.trials.trialTypesMeta(t).randStartPos(i).trials,:),1),2,'gaussian',smoothSpan,'omitnan'),'LineWidth',lineWidth,'Color',color)        
        %plotQuartiles(Xax,binVel(sData.trials.trialTypesMeta(t).randStartPos(i).trials,:),color,smoothSpan)
    end
    axis([corridorStart corridorLength 0 80]);
    xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
    ax = gca;
    ax.TickDir = 'out';
    ylabel('Velocity (cm/s)');
    
    sp2 = subplot(nTrialTypes,2,t*2);
    hold on
    
    rectangle('Position',[corridorStart+1,0.05,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
    rectangle('Position',[rewardZone,0.0025,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
    if stimTrial
    rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
    end
    
    for i = 1:1:4
        color = myMap(i,:);
        plot(Xax,smoothdata(nanmedian(lickFreqInBin(sData.trials.trialTypesMeta(t).randStartPos(i).trials,:),1),2,'gaussian',smoothSpan,'omitnan'),'LineWidth',lineWidth,'Color',color)
        %plotQuartiles(Xax,lickFreqInBin(sData.trials.trialTypesMeta(t).randStartPos(i).trials,:),color,smoothSpan)
    end
    axis([corridorStart corridorLength 0 8]);
    text(120,8.5,['LMI: ' num2str(sData.trials.trialTypesMeta(t).landModInd,2)],'FontSize',10,'FontWeight','bold')
    xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
    ax = gca;
    ax.TickDir = 'out';

    %ylabel('Lick before reward (lick/cm)');
    ylabel('Lick frequency (Hz)');
    
    legend(sp2,plotLegend,'Location','northwest')



end
suptitle([sData.sessionInfo.sessionID newline])


if ~exist([filePath 'landModCurves'], 'dir')
    mkdir([filePath 'landModCurves']);
end

%saveas(gcf,strcat(fullfile(filePath,[sData.sessionInfo.sessionID '_3_landModCurves']),'.png'));
saveas(gcf,fullfile(filePath,'landModCurves',[sData.sessionInfo.sessionID '_landModCurves.png']));

close(gcf)




end
