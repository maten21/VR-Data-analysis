function [] = mapStabilityQC(sData,sDataDir)

fov = 1; 

B1OddEven = sData.imdata.correlations.xEvenOdd(1).maxCorrCoef;
B1B5 = sData.imdata.correlations.xBlocks.maxCorrCoef.B1xB5;

figure('Color','white','Position',[0 0 300 400]);
hold on

Xax = sData.stats.sessionAvs.plotXAxis;
start = find(Xax == 150)+1;
smoothSpan = 5; 


myMap = lines;


plot(Xax,smoothdata(B1OddEven,'gaussian',smoothSpan),'LineWidth',3)
plot(Xax(start:end),smoothdata(B1B5(start:end),'gaussian',smoothSpan),'LineWidth',3)

xlim([0 Xax(end)])
ylim([0 1])

% legend({'Fam. B1 Odd x Even','Fam. B1 x B5'},'Location','east')
text(60,0.35,['Mean = ' num2str(mean(B1OddEven),2)],"FontSize",18,"Color",myMap(1,:))
text(60,0.25,['Mean = ' num2str(mean(B1B5(start:end)),2)],"FontSize",18,"Color",myMap(2,:))


if mean(mean(B1OddEven)) < 0.75
    text(20,0.97,'B1 Odd x Even',"FontSize",13,"Color",myMap(1,:))
    text(20,0.9,'B1 x B5',"FontSize",13,"Color",myMap(2,:))
else
    text(20,0.7,'B1 Odd x Even',"FontSize",13,"Color",myMap(1,:))
    text(20,0.63,'B1 x B5',"FontSize",13,"Color",myMap(2,:))
end

ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');

ylabel('\fontsize{13}Track position (cm)'); 
title('Familiar map stability')

%saveas(gcf,strcat(fullfile(sDataDir,['Fov' num2str(fov) '-' sData.imdata(fov).fovLocation '-famMapStability']),'.png'));
saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '-Fov' num2str(fov) '-' sData.imdata(fov).fovLocation '-famMapStability']),'.png'));

end
