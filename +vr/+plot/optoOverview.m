% Optogenetics light-ON/light-OFF trial figure
%
% optoOverview
%
%
%


% plot dimensions
BCNHeight = (0.5*nBCNTrials/nAllTrials);
NBCHeight = (0.5*nNBCTrials/nAllTrials);
subplotSpacing = 0.05;
relativeWidth = (1 - subplotSpacing*4)/3;

Xax = sessionAvs.plotXAxis;
XaxShort = sessionAvs.plotXAxis(1:(corridorLength-80)/binSize);
binBefRZ = (corridorLength-80)/binSize;

%subplot1:
pos = [subplotSpacing (0.475 + NBCHeight) relativeWidth BCNHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
BCNsubplot = subplot('Position',pos);

imagesc(Xax,1:nBCNTrials,binVel(trialType0,:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,hot);
c.Label.String = 'Speed (virtual cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';
ylabel('Trials');



%Subplot2:
pos = [subplotSpacing 0.4 relativeWidth NBCHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
NBCsubplot = subplot('Position',pos);
imagesc(Xax,1:nNBCTrials,binVel(trialType1,:)) %(1:number of bins;1:number of trials)

c = colorbar;
colormap(gca,hot);
c.Label.String = 'Speed (virtual cm/s)';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 50]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';
ylabel('Trials');


%subplot2.5:

pos = [subplotSpacing 0.125 relativeWidth*0.79 0.2]; %Specify pos as a four-element vector of the form [left bottom width height]. 
VelocitySubplot = subplot('Position',pos);
hold on

rectangle('Position',[-118,0.5,120,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none');
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelBCN,'-',Xax,sessionAvs.avBinVelNBC,'-'); %imagesc(1:BinNumber,1:nNBCTrials,BinVel(TrialType1,:)) %(1:number of bins;1:number of trials)

axis([-120 (binNumber * binSize - 120) 0 60]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');


% for lick plots
myMap = parula(64);
myMap(1,:) = [0 0 0];
%minLick = min(licksInBin)
%maxLick = max(max(licksInBin));
%[0 0 0; parula(5)];
%licksInBin(isnan(licksInBin)) = -1;

%subplot3:
pos = [(subplotSpacing*2+relativeWidth) (0.475 + NBCHeight) relativeWidth BCNHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
LickBCN = subplot('Position',pos);

imagesc(sessionAvs.plotXAxis,1:nBCNTrials,licksInBin(trialType0,:)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';


%Subplot4:
pos = [(subplotSpacing*2+relativeWidth) 0.4 relativeWidth NBCHeight]; %Specify pos as a four-element vector of the form [left bottom width height]. 
LickNBC = subplot('Position',pos);
imagesc(Xax,1:nNBCTrials,licksInBin(trialType1,:)) %(1:number of bins;1:number of trials)

c = colorbar;
colormap(gca,myMap);
c.Label.String = 'Lick count';
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([-0.1 5]); %set limits for color plot, below 1st black, above 2nd white

ax = gca;
ax.TickDir = 'out';






Fig2 = figure('Color','white','Position',[0 0 1200 500]);


VelocitySubplotBCN = subplot(2,3,1);
hold on

rectangle('Position',[-118,0.5,118,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); %R Z
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelBCNLightOff,Xax,sessionAvs.avBinVelBCNLightOn);

axis([-120 (rewardZone+80) 0 60]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');
text((rewardZone-viewDistance),40,'\leftarrow cue visible on horizon');
text(0,45,'\leftarrow corridor starts');
text(-118,50,'"black box"');
text(rewardZone,50,'RZ');
%legend('muscimol','no injection','saline','Location','northwest')

VelocitySubplotNBC = subplot(2,3,4);
hold on

rectangle('Position',[-118,0.5,118,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,20,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(Xax,sessionAvs.avBinVelNBCLightOff,Xax,sessionAvs.avBinVelNBCLightOn);
axis([-120 (rewardZone+80) 0 60]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');


LickBeforeSubplotBCN = subplot(2,3,2);
hold on

rectangle('Position',[-118,0.0025,118,0.6],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,20,0.6],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.0025,0.5,0.6],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(XaxShort,sessionAvs.avLicksInBinBCNLightOff(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinBCNLightOn(1:binBefRZ));
axis([-120 (rewardZone+80) 0 0.6]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick before reward (lick/cm)');
title('Beaconed trials')

LickBeforeSubplotNBC = subplot(2,3,5);
hold on

rectangle('Position',[-118,0.0025,118,0.6],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.0025,20,0.6],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[(rewardZone-viewDistance),0.0025,0.5,0.6],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

plot(XaxShort,sessionAvs.avLicksInBinNBCLightOff(1:binBefRZ),XaxShort,sessionAvs.avLicksInBinNBCLightOn(1:binBefRZ));
axis([-120 (rewardZone+80) 0 0.6]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick before reward (lick/cm)');
title('Non-beaconed trials');


AllLickSubplotBCN = subplot(2,3,3);
semilogy(Xax,sessionAvs.avLicksInBinBCNLightOff,Xax,sessionAvs.avLicksInBinBCNLightOn);
axis([-120 (rewardZone+80) 0.001 10]);
%xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick count (lick/cm)');
legend(['light-off (n=' num2str(nLightOffTrialsBCN) ')'],['light-on (n=' num2str(nLightOnTrialsBCN) ')'],'Location','northwest')


AllLickSubplotNBC = subplot(2,3,6);
semilogy(Xax,sessionAvs.avLicksInBinNBCLightOff,Xax,sessionAvs.avLicksInBinNBCLightOn);
axis([-120 (rewardZone+80) 0.001 10]);
xlabel('Position in unity (virtual cm)');
ax = gca;
ax.TickDir = 'out';
ylabel('Lick count (lick/cm)');
legend(['light-off (n=' num2str(nLightOffTrialsNBC) ')'],['light-on (n=' num2str(nLightOnTrialsNBC) ')'],'Location','northwest')

saveas(gcf,strcat(fullfile(filePath,[sessionID '_2']),'.png'));
close(gcf)

