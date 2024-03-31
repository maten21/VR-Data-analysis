function fig = dffSignals(sData,subset,range)



% plot DFF sinals

%from = min(subset);
%to = max(subset);


if isfield(sData.imdata.meta,'fps')
    fps = sData.imdata.meta.fps;
elseif isfield(sData.imdata.meta,'scanFrameRate') 
    fps = sData.imdata.meta.scanFrameRate;
end

yDff = sData.imdata.roiSignals(2).dff(subset,range);
yDff = okada(okada(yDff,2),2);

yDeconv = sData.imdata.roiSignals(2).deconv(subset,range);
yDeconv(yDeconv>0) = 1;

mymap = lines(numel(subset));

fig = figure('Color','white','Position',[0 0 1200 850]);
hold on
for i = 1:1:numel(subset)

    snr = num2str(sData.imdata.roiStat.signalToNoise(subset(i)),2);
    dffMax = num2str(sData.imdata.roiStat.peakDff(subset(i)),2);
    maxDff = sData.imdata.roiStat.peakDff(subset(i));
    
    x = 1:1:numel(range);
    y = ((yDff(i,:) - min(yDff(i,:)))/maxDff)-i; % ../max(yDff(i,:))-i;
    y2 = yDeconv(i,yDeconv(i,:)>0)*0.8-i;
    
    plot(x,y,'color',mymap(i,:))
    plot(x(yDeconv(i,:)>0), y2,'v','color',mymap(i,:))
    text(double(size(yDff,2)),0.6-i,['      ROI#: ' num2str(subset(i)) newline '      SNR: ' snr newline ' dF/F max: ' dffMax])

end
suptitle([sData.sessionInfo.sessionID '  SNR range: ' num2str(sData.imdata.roiStat.signalToNoise(subset(1)),2) ' - '  num2str(sData.imdata.roiStat.signalToNoise(subset(end)),2)])
%xlim([ ])
rectangle('Position',[numel(range)-60*fps, -10.15, 60*fps, 0.02],'FaceColor',[0 0 0],'EdgeColor','none');
text(numel(range)-30*fps-110,-10.3,'\fontsize{12}1 min')

ylim([-10.5 0.1])
set(gca,'visible','off')






end






