


%% ANALYZE and CLASSIFY ROIs: RELIABLE ACTIVITY CRITERIUM, SMI

% initialize if ftarted with sData
binnedRoisDff = sData.imdata.binnedRoisDff;
binnedRoisDeconv = sData.imdata.binnedRoisDeconv;
binnedRoisDeconvRate = sData.imdata.binnedRoisDeconvRate;
binnedRoisLickAlignedDff = sData.imdata.binnedRoisLickAlignedDff;
binnedRoisLickAlignedDeconv = sData.imdata.binnedRoisLickAlignedDeconv;
binnedRoisLickAlignedDeconvRate = sData.imdata.binnedRoisLickAlignedDeconvRate;

avBinnedRoisDff = sData.imdata.avBinnedRois.avBinnedRoisDff;
avBinnedRoisDeconv = sData.imdata.avBinnedRois.avBinnedRoisDeconv;
avBinnedRoisDeconvRate = sData.imdata.avBinnedRois.avBinnedRoisDeconvRate;
avBinnedRoisLickAlignedDff = sData.imdata.avBinnedRois.avBinnedRoisLickAlignedDff;
avBinnedRoisLickAlignedDeconv = sData.imdata.avBinnedRois.avBinnedRoisLickAlignedDeconv;
avBinnedRoisLickAlignedDeconvRate = sData.imdata.avBinnedRois.avBinnedRoisLickAlignedDeconvRate;


nROIs = sData.imdata.nROIs;
binnedRoisDff = sData.imdata.binnedRoisDeconv;
binnedRoisDeconv = sData.imdata.binnedRoisDeconv;
binnedRoisDeconvRate = sData.imdata.binnedRoisDeconvRate;
allTrials = sData.imdata.nAllTrials;

