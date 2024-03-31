
function sData =  transposeRoiSignals(sData) 

if size(sData.imdata.roiSignals(2).dff,1) ~= sData.imdata.nROIs && size(sData.imdata.roiSignals(2).dff,1) == sData.imdata.nSamples

sData.imdata.roiSignals(2).roif = sData.imdata.roiSignals(2).roif';
sData.imdata.roiSignals(2).npilf = sData.imdata.roiSignals(2).npilf';
sData.imdata.roiSignals(2).dff = sData.imdata.roiSignals(2).dff';
sData.imdata.roiSignals(2).deconv = sData.imdata.roiSignals(2).deconv';
sData.imdata.roiSignals(2).denoised = sData.imdata.roiSignals(2).denoised';
sData.imdata.roiSignals(2).actRateDeconv = sData.imdata.roiSignals(2).actRateDeconv';
end
end




