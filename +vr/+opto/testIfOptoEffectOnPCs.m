

for roi = 1:1:sData.imdata.nROIs

bootStrapData = vr.isSignificantBootstrap(sData.imdata.binned.RoidFF{1, roi}  ,sData.behavior.opto.OptoOnTrialsIndices);

roiMeta(roi,1).p = bootStrapData.p;
roiMeta(roi,1).isSign = bootStrapData.p < 0.05;
%roiMeta(roi,1).bootStrapData;
if floor(roi/10) == ceil(roi/10)
    
    disp(roi);
   
end
end



figure



sum(bootStrapData.amplOpto < bootStrapData.ampl2)/length(bootStrapData.ampl2);
sum(bootStrapData.amplCtrl < bootStrapData.ampl1)/length(bootStrapData.ampl1)


figure
hold on

plot(bootStrapData.ampl1,bootStrapData.ampl2,'.')
plot(bootStrapData.amplCtrl,bootStrapData.amplOpto,'o')

range1 = min([min([bootStrapData.ampl1, bootStrapData.ampl2]), bootStrapData.amplCtrl, bootStrapData.amplOpto])*0.9; 
range2 = max([max([bootStrapData.ampl1, bootStrapData.ampl2]), bootStrapData.amplCtrl, bootStrapData.amplOpto])*1.1; 
min([bootStrapData.ampl1, bootStrapData.ampl2])];
xlim([range1 range2])
ylim([range1 range2])

xlabel('Control ampl.')
ylabel('Opto ampl.')


figure
hold on

histogram(bootStrapData.ampl2)
plot(bootStrapData.amplOpto,500,'o')


figure
hold on

histogram(bootStrapData.ampl1)
plot(bootStrapData.amplCtrl,500,'o')


effectRand = bootStrapData.ampl2-bootStrapData.ampl1;
effect = bootStrapData.amplOpto - bootStrapData.amplCtrl;

sum(effect < effectRand)/length(effectRand)

figure
hold on

histogram(effectRand)
plot(effect,500,'o')


