
% NOT FINISHED, NEEDS A LOT OF WORK



function sData = classifyROIs(sData,sDataDir)


nFOVs =  length(sData.imdata);

for fov = 1:1:nFOVs
    
    sData = classifyROIsFOV(sData, fov, sDataDir);
end

end


function sData = classifyROIsFOV(sData, fov, sDataDir)

%% Create roiMeta structure

nROIs = sData.imdata(fov).nROIs;

if ~isfield(sData.imdata(fov),'roiMeta')
    sData.imdata(fov).roiMeta(1:nROIs) = struct();
end

% SNR, baseline...
% sData.imdata(i).roiMeta(roi).SNR = ;

for roi = 1:1:nROIs
    sData.imdata(fov).roiMeta(roi).SNR = sData.imdata(fov).roiStat.signalToNoise(roi);
    sData.imdata(fov).roiMeta(roi).peakDff = sData.imdata(fov).roiStat.peakDff(roi);
    sData.imdata(fov).roiMeta(roi).activityLevel = sData.imdata(fov).roiStat.activityLevel(roi);
end




nTrialTypes = length(sData.trials.contextsMeta);

%nTrialTypes = 4;    
%trialType = 1; % set trial or contect type


for trialType = 1:1:nTrialTypes
posTunCurves = sData.imdata(fov).avBinnedRois.avBinnedRoisDff{trialType};
binnedRoisDff = sData.imdata(fov).binnedRoisDff;
trials = sData.trials.contextsMeta(trialType).trials;
% nROIs = sData.imdata(fov).nROIs;
try 
binSize = sData.imdata(fov).binSize;
catch
binSize = sData.behavior.trialMatrices.meta.binSize; 
end

try
binNumber = sData.imdata(fov).binNumber;    
catch
    binNumber = sData.behavior.trialMatrices.meta.binNumber;
end



%allTrials = sData.imdata(fov).nAllTrials;



smoothSpan = 10;
firstThreshold = 0.3;  % detect activities in the position tuning curves
secondThreshold = 2.5;
fieldWidthThreshold = 12;   % cm should be a multiplicate of the binSize
reliabilityThreshold = 0.5; % reliability among trials 



for roi = 1:1:nROIs
    roiMetaArray(roi,trialType).placeCell = true;
end


%% --- #0: quality check --- %%
% exclude data with poor SNR 
 

%% --- #1: identify potential fields  --- %%

% Initial smoothing of tuning curves
posTunCurves =  smoothdata(posTunCurves,2,'gaussian',smoothSpan);

% determine potential fields: min activity > first threshold
[peakActivity, peakActPos] = max(posTunCurves'); % [maxVal, index] = max(A)

% peak align data based on position tuning curve (use same for all trials)
trackCenter = round(binNumber/2);
shiftPeak = trackCenter - peakActPos;
peakAlignedPosTunCurves = posTunCurves;
for roi = 1:1:nROIs
    shiftArray = 1:1:binNumber;    
    if shiftPeak(roi) > 0
        shiftArray = [shiftArray(end-shiftPeak(roi)+1:end), shiftArray(1:end-shiftPeak(roi))];
    elseif shiftPeak(roi) < 0
        shiftArray = [shiftArray(abs(shiftPeak(roi))+1:end), shiftArray(1:abs(shiftPeak(roi)))];
    end
    peakAlignedPosTunCurves(roi,:) = posTunCurves(roi,shiftArray);
end
posTunCurves = peakAlignedPosTunCurves;




baseline = quantile(posTunCurves,0.1,2);
potFields = posTunCurves;

for roi = 1:1:nROIs
    
    potFields(roi,:) = (posTunCurves(roi,:) - baseline(roi)) > (peakActivity(roi)-baseline(roi)) * firstThreshold;
    snr(roi) = sData.imdata(fov).roiMeta(roi).SNR;
end

% Find the start and end of potential fields, than generate roi matrix where
% the field length is written in the starting bin 
diffRois = zeros(nROIs,binNumber);
diffRois(:,2:binNumber) = diff(potFields,1,2);

fieldStarts = zeros(nROIs,binNumber); % Will contain the field lengths at the starting bin
for roi = 1:1:nROIs
    tempStarts = find(diffRois(roi,:) == 1);
    tempEnds = find(diffRois(roi,:) == -1);
    
    if numel(tempStarts) == numel(tempEnds)
        fieldStarts(roi,tempStarts) = tempEnds - tempStarts;
        
    elseif numel(tempStarts) > numel(tempEnds) % corrects for field at the end of the track
        tempEnds = [tempEnds binNumber];
        fieldStarts(roi,tempStarts) = tempEnds - tempStarts;
        
    else                                       % corrects for field at the beginning of the track
        tempStarts = [1 tempStarts];
        fieldStarts(roi,tempStarts) = tempEnds - tempStarts;
    end
    
    clear('tempStarts','tempEnds');
end

% set to zero those potential fields that are below the threshold
fieldStarts(fieldStarts < fieldWidthThreshold/binSize) = 0;

for roi = 1:1:nROIs
    if max(fieldStarts(roi,:)) == 0
        roiMetaArray(roi,trialType).placeCell = false;
    end
    roiMetaArray(roi,trialType).nFields = numel(find(fieldStarts(roi,:)));
end


%%% --- #2: mean in field activity must be 3x lagher than the out of field activity  --- %%%

% generate matrix containing 1 at in fiels 0 at out of fields
fields = fieldStarts;

for roi = 1:1:nROIs
    if roiMetaArray(roi,trialType).placeCell == true
        tempRoi = fieldStarts(roi,:);
        tempStarts = find(tempRoi);
        tempEnds = tempStarts + tempRoi(find(tempRoi))-1;
        for f = 1:1:roiMetaArray(roi,trialType).nFields
            fields(roi,tempStarts(f):tempEnds(f)) = 1;
        end
        clear('tempRoi','tempStarts','tempEnds');
    end
end

% exclude cells that do not meet the criterion
counter = 0;
for roi = 1:1:nROIs
    if roiMetaArray(roi,trialType).placeCell == true
        
        tempFields = fields(roi,:);
        tempRoi = posTunCurves(roi,:);
        
        if nanmean(tempRoi(find(tempFields))) < nanmean(tempRoi(find(tempFields == 0)))*secondThreshold
            roiMetaArray(roi,trialType).placeCell = false;
            counter = counter + 1; 
        end
        clear('tempRoi','tempFields');
    end
end


%%% --- #3: Reliability: place cells must be active in 30% of trials --- %%%
% activity peak should be in field



% This is how to shift the entire trial matrix in circular manner
%{
testMatr = testMatr0';
A = testMatr(:)';
A2 = [A(20:end), A(1:19)];
testMatr2 = reshape(A2',size(testMatr))';
%}




% find the peak activity position for each ROIs in all trials
peakActPosInTrials = zeros(nROIs,numel(trials));
binnedRoisDffPeakAligned = [];
for roi = 1:1:nROIs
    
    M = binnedRoisDff(:,:,roi)'; A = M(:)';
    if shiftPeak(roi) < 0
        A2 = [A(abs(shiftPeak(roi))+1:end), A(1:abs(shiftPeak(roi)))];
    elseif shiftPeak(roi) > 0
        A2 = [A(end-shiftPeak(roi)+1:end), A(1:end-shiftPeak(roi))];
    else
        A2 = A;
    end
    shiftedMatr = reshape(A2',size(M))';
    binnedRoisDffPeakAligned(:,:,roi) = shiftedMatr;
    for t = 1:1:numel(trials)
        if numel(find(binnedRoisDff(trials(t),:,roi) == max(binnedRoisDff(trials(t),:,roi)))) == 1 %% only consider if there is only one solution, i.e. exclude inactive trials with only zeros
            maxValLocation = find(shiftedMatr(trials(t),:) == max(shiftedMatr(trials(t),:)));
            peakActPosInTrials(roi,t) = maxValLocation(1);
            %if numel(find(binnedRoisDff(trials(t),:,roi) == max(binnedRoisDff(trials(t),:,roi)))) == 1 %% only consider if there is only one solution, i.e. exclude inactive trials with only zeros
            %peakActPosInTrials(roi,t) = find(binnedRoisDff(trials(t),:,roi) == max(binnedRoisDff(trials(t),:,roi)));
        end
    end
end


% Check if the peak activity is member of the in field array
counter = 0;
for roi = 1:1:nROIs
    if roiMetaArray(roi,trialType).placeCell == true
        tempInField = find(fields(roi,:));
        if sum(ismember(peakActPosInTrials(roi,:),tempInField)) < numel(trials)*reliabilityThreshold
            roiMetaArray(roi,trialType).placeCell = false;
            counter = counter + 1;
        end
        clear('tempInField');
    end
end

% Clear false place fields
for roi = 1:1:nROIs
    if roiMetaArray(roi,trialType).placeCell == false
        roiMetaArray(roi,trialType).nFields = 0;
    end
end

% Extract and count place cell ROI indexes
placeCells = [];
for roi = 1:1:nROIs
    if roiMetaArray(roi,trialType).placeCell
        placeCells = [placeCells roi];     
    end
    
    sData.imdata(fov).roiMeta(roi).placeCell(trialType) = roiMetaArray(roi,trialType).placeCell;
    sData.imdata(fov).roiMeta(roi).nFields(trialType) = roiMetaArray(roi,trialType).nFields;
end
nPlaceCells = numel(placeCells);
placeCellFraction = nPlaceCells/nROIs;


%sData.trials.contextsMeta(trialType).placeCells(fov) = placeCells;
%sData.trials.contextsMeta(trialType).nPlaceCells(fov) = nPlaceCells;
%sData.trials.contextsMeta(trialType).placeCellFraction(fov) = placeCellFraction;

blockData(trialType).placeCells = placeCells;
blockData(trialType).nPlaceCells = nPlaceCells;
blockData(trialType).placeCellFraction = placeCellFraction;






%% Control Figure

Xax = sData.stats.sessionAvs(1).plotXAxis;
cMin = quantile(posTunCurves(:),0.01);
cMax = quantile(posTunCurves(:),0.99);
cLabel = 'Dff';

trialsOdd = sData.trials.contextsMeta(trialType).trials(1:2:sData.trials.contextsMeta(trialType).nTrials);
trialsEven = sData.trials.contextsMeta(trialType).trials(2:2:sData.trials.contextsMeta(trialType).nTrials);

posTunCurvesOdd =  nanmean(sData.imdata(fov).binnedRoisDff(trialsOdd,:,:)); 
posTunCurvesOdd = permute(posTunCurvesOdd,[3 2 1]);

posTunCurvesEven =  nanmean(sData.imdata(fov).binnedRoisDff(trialsEven,:,:)); 
posTunCurvesEven = permute(posTunCurvesEven,[3 2 1]);
posTunCurvesEven =  smoothdata(posTunCurvesEven,2,'gaussian',5);

%posTunCurves = sData.imdata(fov).avBinnedRois.avBinnedRoisDff{trialType};
%posTunCurves =  smoothdata(posTunCurves,2,'gaussian',smoothSpan);


figure('Color','white','Position',[0 0 800 400])

sorted = vr.sortROIs(posTunCurvesOdd,placeCells,5);

subplot(1,2,1)
imagesc(Xax,1:numel(sorted),posTunCurvesEven(sorted,:))
%title('Place cells') 
t = title('Place cells');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
c.Label.String = cLabel;
caxis([cMin cMax]);


sorted = vr.sortROIs(posTunCurvesOdd,setdiff(1:nROIs,placeCells),5);

subplot(1,2,2)
imagesc(Xax,1:numel(sorted),posTunCurvesEven(sorted,:))
%title('Non-place cells') 
t = title('Non-place cells');
t.FontWeight = 'normal';
t.FontSize = 14;
c = colorbar;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
c.Label.String = cLabel;
caxis([cMin cMax]);

suptitle([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' - ' sData.trials.vrContextProtocols(trialType).name num2str(trialType)])

saveas(gcf,strcat(fullfile(sDataDir,[sData.sessionInfo.sessionID(1:17) '_PlaceCellsSorted_dff_' 'Fov' num2str(fov) '_' sData.imdata(fov).fovLocation '_Block' num2str(trialType)]),'.png'));



%% Find inducation lap


% sorted = vr.sortROIs(posTunCurves,placeCells,5);
% posTunCurves(sorted,:)


% Calculate mean activity within field for each trials

peakPos = nan(1,nPlaceCells);

data = zeros(nPlaceCells,numel(trials));
windowWidth = 2; % later probably a better definition is needed maybe concider entire place field % peakPos - window : peakPos : peakPos + window)
for i = 1:1:nPlaceCells

    roi = placeCells(i);
    
    peakPos(i) = find(max(posTunCurves(roi,:)') == posTunCurves(roi,:)');
    
    if floor(peakPos(i)) - windowWidth < 1
        avWindow = 1:1:(floor(peakPos(i)) + windowWidth+1);
    elseif floor(peakPos(i)) + windowWidth > binNumber
        avWindow = (floor(peakPos(i)) - windowWidth):1:binNumber;
    else
        avWindow = (floor(peakPos(i)) - windowWidth):1:(floor(peakPos(i)) + windowWidth);
    end
    
       
    % data(i,:) = nanmean(sData.imdata(fov).binnedRoisDff(trials,avWindow,roi),2);
    data(i,:) = nanmean(binnedRoisDffPeakAligned(trials,avWindow,roi),2); % peaks are aligned to the center bin68
end

%   figure; plot(nanmean(sData.imdata(fov).binnedRoisDff(trials,avWindow,roi),2));
    
    
    
 
 potFieldStarts = data;

 peakActivity = nanmean(data,2); % sorted place cells only
 baseline = quantile(posTunCurves,0.1,2);
 baseline = baseline(placeCells);  % keep the otiginal out of field

 for roi = 1:1:nPlaceCells
    
    potFieldStarts(roi,:) = (data(roi,:) - baseline(roi)) > (peakActivity(roi)-baseline(roi)) * firstThreshold;
    
end
 



    
% Find the start and end of potential fields, than generate roi matrix where the field length is written in the starting bin 
intWindow = 5;

diffRois = zeros(nPlaceCells,numel(trials)+2); 
diffRois = diff([false(nPlaceCells,1) potFieldStarts false(nPlaceCells,1)],1,2); %zeros are added to the beginning and end to see place fields start in the first trial

inductionTrial = zeros(nPlaceCells,1); 
for roi = 1:1:nPlaceCells
    tempStarts = find(diffRois(roi,:) == 1);
    tempEnds = find(diffRois(roi,:) == -1);

    
    for i = 1:1:numel(tempStarts)
        
        if tempStarts(i)+intWindow-1 > size(potFieldStarts,2)
            if sum(potFieldStarts(roi,tempStarts(i):size(potFieldStarts,2))) > 1
                inductionTrial(roi,trialType) = tempStarts(i);
                break
            else
                inductionTrial(roi,trialType) = 0;
                break
            end
        elseif sum(potFieldStarts(roi,tempStarts(i):tempStarts(i)+intWindow-1)) > 1
            
            inductionTrial(roi,trialType) = tempStarts(i);
            break
        end
        
    end
        
        
        %clear('tempStarts','tempEnds');
end


blockData(trialType).inductionTrial = inductionTrial(:,trialType);
blockData(trialType).peakPos = peakActPos(placeCells)';



end

sData.imdata(fov).placeCells = blockData;



%% ACTIVE TRIAL FRACTION


peakFreqs(1:nROIs) = nan;
activeTrialFraction = [];

for roi = 1:1:nROIs   % ACTIVE TRIAL FRACTION
    %activeTrials              = sum(nanmean(sData.imdata(fov).binnedRoisDeconv(:,:,roi)') > 0) / sData.imdata(fov).nAllTrials;
    sData.imdata(fov).roiMeta(roi).activeTrialFraction = sum(nanmean(sData.imdata(fov).binnedRoisDeconv(:,:,roi)') > 0) / size(sData.imdata(fov).binnedRoisDff,1);
    activeTrialFraction(roi) = sum(nanmean(sData.imdata(fov).binnedRoisDeconv(:,:,roi)') > 0) / size(sData.imdata(fov).binnedRoisDff,1);            
    %data = binnedRoisDeconvRate(:,:,roi);
    %        roiMeta(roi).peakFreq = quantile(data(:),[0.99]);
    %        peakFreqs(roi) = quantile(data(:),[0.99]);
    
   for t = 1:1:length(sData.trials.contextsMeta)
    trials = sData.trials.contextsMeta(t).trials;
        sData.imdata(fov).roiMeta(roi).activeTrialFractionInBlock(t) = sum(nanmean(sData.imdata(fov).binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
       % activeTrialFractionInBlock(roi,t)                           = sum(nanmean(sData.imdata(fov).binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
       %      roiMeta(roi).activeTrialFractionInBlock(t) = sum(nanmean(sData.imdata(fov).binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
        
     %       data = binnedRoisDeconvRate(trials,:,roi);
     %       roiMeta(roi).peakFreqInBlock(t) = quantile(data(:),[0.99]);
            % peakFreqsCtrl(roi) = quantile(data(:),[0.99]);
    
            % roiMeta(roi).controlSMI = vr.getSMI(binnedRoisDeconvRate(trials,31:165,roi)); % Calculate SMI
   end
   
end

% Active ROIs
threshold = 0.3;
%freqThreshold = 1.5;

% activeROIs = intersect(find(activeTrialFractionCtrl > threshold),find(peakFreqsCtrl > freqThreshold));    
sData.imdata(fov).activeROIs = find(activeTrialFraction > threshold); 
sData.imdata(fov).inactiveROIs = setdiff(1:nROIs,sData.imdata(fov).activeROIs)';

close all
end


% Earlier version, I am not sure if it has ever worked. Keep temporarily!
%{
function sData = classifyROIsFOV(sData, i)



%% Calculate metadata for rois (important for pre selection, quality control)

nROIs = sData.imdata(i).nROIs;

if ~isfield(sData.imdata(i),'roiMeta')
    sData.imdata(i).roiMeta(1:nROIs) = struct();
end

% SNR, baseline...
% sData.imdata(i).roiMeta(roi).SNR = ;

for roi = 1:1:nROIs
    sData.imdata(i).roiMeta(roi).SNR = sData.imdata(i).roiStat.signalToNoise(roi);
    sData.imdata(i).roiMeta(roi).peakDff = sData.imdata(i).roiStat.peakDff(roi);
    sData.imdata(i).roiMeta(roi).activityLevel = sData.imdata(i).roiStat.activityLevel(roi);
end




% initialize if ftarted with sData
binnedRoisDff = sData.imdata(i).binnedRoisDff;
binnedRoisDeconv = sData.imdata(i).binnedRoisDeconv;
binnedRoisDeconvRate = sData.imdata(i).binnedRoisDeconvRate;
%binnedRoisLickAlignedDff = sData.imdata(i).binnedRoisLickAlignedDff;
%binnedRoisLickAlignedDeconv = sData.imdata(i).binnedRoisLickAlignedDeconv;
%binnedRoisLickAlignedDeconvRate = sData.imdata(i).binnedRoisLickAlignedDeconvRate;

avBinnedRoisDff = sData.imdata(i).avBinnedRois.avBinnedRoisDff;
avBinnedRoisDeconv = sData.imdata(i).avBinnedRois.avBinnedRoisDeconv;
avBinnedRoisDeconvRate = sData.imdata(i).avBinnedRois.avBinnedRoisDeconvRate;
%avBinnedRoisLickAlignedDff = sData.imdata(i).avBinnedRois.avBinnedRoisLickAlignedDff;
%avBinnedRoisLickAlignedDeconv = sData.imdata(i).avBinnedRois.avBinnedRoisLickAlignedDeconv;
%avBinnedRoisLickAlignedDeconvRate = sData.imdata(i).avBinnedRois.avBinnedRoisLickAlignedDeconvRate;



%binnedRoisDff = sData.imdata(i).binnedRoisDff;
%binnedRoisDeconv = sData.imdata(i).binnedRoisDeconv;
%binnedRoisDeconvRate = sData.imdata(i).binnedRoisDeconvRate;




%% Test place cell criteria (based on Mao's criteria)
% store boolean values in "roiMeta(roi).placeCell" about wether a roi passes place cell criteria or not
%
% CRITERIA FOR PLACE CELLS
% 1) Dombeck 2010: Potential place fields were first identified as contiguous regions of
%    this plot in which all of the points were greater than 25% of the difference between
%    the peak ?F/F value (for all 80 bins) and the baseline value (mean of the lowest
%    20 out of 80 ?F/F values). I Set a 30% threshold as Mao.
% 2) Dombeck 2010: The mean in-field ?F/F value must be more than three times the mean out-of-field ?F/F value.
% 3) Significant calcium transients must be present >30% of the time the mouse spent in the place field.
%
%    The mean in-field activity of the position activity map for the
%    deconvolved signal must be at least 3 times larger than the mean
%    out-of-field activity.
% 3) The peak firing for an individual lap has to be within the field as
%    classified with the position activity map for at least one-third of all
%    laps.
% 4) UNDER DEVELOPMENT The spatial information (SI) ....
%
%
% SPECIFY PARAMETERS
% SPECIFY WHICH TRIAL SUBSET AND POSITION TUNING CURVES SHOULD BE USED:


posTunCurves = avBinnedRoisDeconv{1};

%{
% Correct for non random reward start
for roi = 1:1:nROIs

    posTunCurves(roi,[1:30 size(posTunCurves,2)-29:size(posTunCurves,2)]) = nanmean(posTunCurves(roi,31:(size(posTunCurves,2)-30)));
end
%}

trialSubset = sData.trials.trialTypesMeta(1).trials;
%meta.posTunCurvesUsed = 'avBinnedRoisBCNLightOffDeconv'; %UPDATE MANUALLY!!!
%meta.trialSubsetUsed = 'lightOffBCNTrials'; %UPDATE MANUALLY!!!


smoothSpan = 10;
firstThreshold = 0.3;       % detect activities in the position tuning curves
fieldWidthThreshold = 12;   % cm should be a multiplicate of the binSize
reliabilityThreshold = 0.3; % reliability among trials

binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
allTrials = sData.imdata(i).nAllTrials;


for roi = 1:1:nROIs
    sData.imdata(i).roiMeta(roi).placeCell = true;
end

%%% --- #0: quality check --- %%%
% exclude data with poor SNR


%%% --- #1: identify potential fields  --- %%%

% Initial smoothing of tuning curves
posTunCurves =  smoothdata(posTunCurves,2,'gaussian',smoothSpan);

% determine potential fields: min activity > first threshold
[peakActivity, peakActPos] = max(posTunCurves'); % [maxVal, index] = max(A)
potFields = posTunCurves;

for roi = 1:1:nROIs
    
    potFields(roi,:) = posTunCurves(roi,:) > peakActivity(roi) * firstThreshold;
    
end

% Find the start and end of potential fields, than generate roi matrix where
% the field length is written in the starting bin
diffRois = zeros(nROIs,binNumber);
diffRois(:,2:binNumber) = diff(potFields,1,2);

fieldStarts(nROIs,binNumber) = zeros; % Will contain the field lengths at the starting bin
for roi = 1:1:nROIs
    tempStarts = find(diffRois(roi,:) == 1);
    tempEnds = find(diffRois(roi,:) == -1);
    
    if numel(tempStarts) == numel(tempEnds)
        fieldStarts(roi,tempStarts) = tempEnds - tempStarts;
        
    elseif numel(tempStarts) > numel(tempEnds) % corrects for field at the end of the track
        tempEnds = [tempEnds binNumber];
        fieldStarts(roi,tempStarts) = tempEnds - tempStarts;
        
    else                                       % corrects for field at the beginning of the track
        tempStarts = [1 tempStarts];
        fieldStarts(roi,tempStarts) = tempEnds - tempStarts;
    end
    
    clear('tempStarts','tempEnds');
end

% set to zero those potential fields that are below the threshold
fieldStarts(fieldStarts < fieldWidthThreshold/binSize) = 0;

for roi = 1:1:nROIs
    if max(fieldStarts(roi,:)) == 0
        sData.imdata(i).roiMeta(roi).placeCell = false;
    end
    sData.imdata(i).roiMeta(roi).nFields = numel(find(fieldStarts(roi,:)));
end


%%% --- #2: mean in field activity must be 3x lagher than the out of field activity  --- %%%

% generate matrix containing 1 at in fiels 0 at out of fields
fields = fieldStarts;

for roi = 1:1:nROIs
    if sData.imdata(i).roiMeta(roi).placeCell == true
        tempRoi = fieldStarts(roi,:);
        tempStarts = find(tempRoi);
        tempEnds = tempStarts + tempRoi(find(tempRoi))-1;
        for f = 1:1:sData.imdata(i).roiMeta(roi).nFields
            fields(roi,tempStarts(f):tempEnds(f)) = 1;
        end
        clear('tempRoi','tempStarts','tempEnds');
    end
end

% exclude place cells that do not meet the criterion
counter = 0;
for roi = 1:1:nROIs
    if sData.imdata(i).roiMeta(roi).placeCell == true
        
        tempFields = fields(roi,:);
        tempRoi = posTunCurves(roi,:);
        
        if nanmean(tempRoi(find(tempFields))) < nanmean(tempRoi(find(tempFields == 0)))*3
            sData.imdata(i).roiMeta(roi).placeCell = false;
            counter = counter + 1;
        end
        clear('tempRoi','tempFields');
    end
end


%%% --- #3: Reliability: place cells must be active in 30% of trials --- %%%
% activity peak should be in field

% find the peak activity position for each ROIs in all trials
peakActPosInTrials(nROIs,allTrials) = 0;

for roi = 1:1:nROIs
    for t = 1:1:allTrials
        if numel(find(binnedRoisDeconv(t,:,roi) == max(binnedRoisDeconv(t,:,roi)))) == 1 %% only consider if there is only one solution, i.e. exclude inactive trials with only zeros
            peakActPosInTrials(roi,t) = find(binnedRoisDeconv(t,:,roi) == max(binnedRoisDeconv(t,:,roi)));
        end
    end
end


% Check if the peak activity is member of the in field array
counter = 0;
for roi = 1:1:nROIs
    if sData.imdata(i).roiMeta(roi).placeCell == true
        tempInField = find(fields(roi,:));
        if sum(ismember(peakActPosInTrials(roi,trialSubset),tempInField)) < numel(trialSubset)*reliabilityThreshold
            sData.imdata(i).roiMeta(roi).placeCell = false;
            counter = counter + 1;
        end
        clear('tempInField');
    end
end

% Clear false place fields
for roi = 1:1:nROIs
    if sData.imdata(i).roiMeta(roi).placeCell == false
        sData.imdata(i).roiMeta(roi).nFields = 0;
    end
end

% Extract and count place cell ROI indexes
placeCells = [];
for roi = 1:1:nROIs
    if sData.imdata(i).roiMeta(roi).placeCell
        placeCells = [placeCells roi];
    end
end
nPlaceCells = numel(placeCells);
placeCellFraction = nPlaceCells/nROIs;

%sData.imdata(i).roiMeta = roiMeta;




%% ANALYZE and CLASSIFY ROIs: RELIABLE ACTIVITY CRITERIUM, SMI


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Xax = ([1:size(binnedRoisLickAlignedDeconvRate,2)]-161)*binSize;

peakFreqs(1:nROIs) = nan;

for roi = 1:1:nROIs   % ACTIVE TRIAL FRACTION
    sData.imdata(i).roiMeta(roi).activeTrials = sum(nanmean(binnedRoisDeconv(:,:,roi)') > 0) / allTrials;
    activeTrialFraction(roi,1) = sum(nanmean(binnedRoisDeconv(:,:,roi)') > 0) / allTrials;
    
    data = binnedRoisDeconvRate(:,:,roi);
    sData.imdata(i).roiMeta(roi).peakFreq = quantile(data(:),[0.99]);
    peakFreqs(roi) = quantile(data(:),[0.99]);
    
    for t = 1:1:numel(sData.trials.contextsMeta)
        trials = sData.trials.contextsMeta(t).trials;
        sData.imdata(i).roiMeta(roi).activeTrialsInBlock(t) = sum(nanmean(binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
        activeTrialFractionInBlock(roi,1) = sum(nanmean(binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
        sData.imdata(i).roiMeta(roi).activeTrialFractionInBlock(t) = sum(nanmean(binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
        
        data = binnedRoisDeconvRate(trials,:,roi);
        sData.imdata(i).roiMeta(roi).peakFreqInBlock(t) = quantile(data(:),[0.99]);
        % peakFreqsCtrl(roi) = quantile(data(:),[0.99]);
        
        % sData.imdata(i).roiMeta(roi).controlSMI = vr.getSMI(binnedRoisDeconvRate(trials,31:165,roi)); % Calculate SMI
    end
    %{
if size(sData.trials.trialTypesMeta,2) > 1
    trials = sData.trials.trialTypesMeta(2).trials;
        sData.imdata(i).roiMeta(roi).activeTrialsHB = sum(nanmean(binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
        activeTrialFractionHB(roi,1) = sum(nanmean(binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
    
            data = binnedRoisDeconvRate(trials,:,roi);
            peakFreqsHB(roi) = quantile(data(:),[0.99]);
    
    trials = sData.trials.trialTypesMeta(3).trials;
        sData.imdata(i).roiMeta(roi).activeTrialsOTP = sum(nanmean(binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);
        activeTrialFractionOTP(roi,1) = sum(nanmean(binnedRoisDeconv(trials,:,roi)') > 0) / numel(trials);

            data = binnedRoisDeconvRate(trials,:,roi);
            peakFreqsOTP(roi) = quantile(data(:),[0.99]);
end
    %}
    
end

for roi = 1:1:nROIs  % Shuffling
    for t = 1:1:numel(sData.trials.contextsMeta)
        
        trials = sData.trials.contextsMeta(t).trials;
        
        [~, randPeakAmpls]  = vr.getSMI(binnedRoisDff(trials,:,roi)); % Calculate SMI
        
        sData.imdata(i).roiMeta(roi).shuffleMeanDff(t) = mean(randPeakAmpls);
        sData.imdata(i).roiMeta(roi).shuffleStdDff(t) = std(randPeakAmpls);
        
        [~, randPeakAmpls] = vr.getSMI(binnedRoisDeconvRate(trials,:,roi)); % Calculate SMI
        
        sData.imdata(i).roiMeta(roi).shuffleMean(t) = mean(randPeakAmpls);
        sData.imdata(i).roiMeta(roi).shuffleStd(t) = std(randPeakAmpls);
        
    end
end

% Save Track aligned and lick aligned peaks of the position tuning curves,
% SMI Z-score
smoothSpan = 9;
lickCoupledROIs = NaN;
lickCoupledROIsDff = NaN;
i = 1;
j = 1;
for roi = 1:1:nROIs % Calculate SMI for both lick- and track- aligned
    
    for t = 1:1:numel(sData.trials.contextsMeta)
        % Deconv Track
        tunCurve = smoothdata(avBinnedRoisDeconvRate{1,t}(roi,:),2,'gaussian',smoothSpan);
        sData.imdata(i).roiMeta(roi).peakDeconv(t) = (max(tunCurve) - min(tunCurve)); % / mean(tunCurve);
        peakBin = find(tunCurve == max(tunCurve));
        sData.imdata(i).roiMeta(roi).peakBinDeconv(t) = peakBin(1);
        peakBinDeconv(roi,t) = peakBin(1);
        
        % Dff Track
        tunCurve = smoothdata(avBinnedRoisDff{1,t}(roi,:),2,'gaussian',smoothSpan);
        tunCurve(isnan(tunCurve)) = 0;
        sData.imdata(i).roiMeta(roi).peakDff(t) = (max(tunCurve) - min(tunCurve)); % / mean(tunCurve);
        peakBin = find(tunCurve == max(tunCurve));
        sData.imdata(i).roiMeta(roi).peakBinDff(t) = peakBin(1);
        peakBinDff(roi,t) = peakBin(1);
        
        
        % Deconv
        sData.imdata(i).roiMeta(roi).SMIDeconv(t) = abs(roiMeta(roi).shuffleMean(t) - sData.imdata(i).roiMeta(roi).peakDeconv(t)) / sData.imdata(i).roiMeta(roi).shuffleStd(t);
        SMIDeconv(roi,t) = sData.imdata(i).roiMeta(roi).SMIDeconv(t);
        %roiMeta(roi).lickAlCtrlSMI = abs(roiMeta(roi).shuffleMean - sData.imdata(i).roiMeta(roi).lickAlignedPeak) / sData.imdata(i).roiMeta(roi).shuffleStd;
        %    lickAlCtrlSMI(roi,1) = sData.imdata(i).roiMeta(roi).lickAlCtrlSMI;
        
        % Dff
        sData.imdata(i).roiMeta(roi).SMIDff(t) = abs(roiMeta(roi).shuffleMeanDff(t) - sData.imdata(i).roiMeta(roi).peakDff(t)) / sData.imdata(i).roiMeta(roi).shuffleStdDff(t);
        SMIDff(roi,t) = sData.imdata(i).roiMeta(roi).SMIDff(t);
        %roiMeta(roi).lickAlCtrlSMIDff = abs(roiMeta(roi).shuffleMeanDff - sData.imdata(i).roiMeta(roi).lickAlignedPeakDff) / sData.imdata(i).roiMeta(roi).shuffleStdDff;
        %   lickAlCtrlSMIDff(roi,1) = sData.imdata(i).roiMeta(roi).lickAlCtrlSMIDff;
    end
    
    %{
    if  sData.imdata(i).roiMeta(roi).activeTrialFractionInBlock > 0.3 && sData.imdata(i).roiMeta(roi).lickAlCtrlSMI > 3  && sData.imdata(i).roiMeta(roi).lickAlignedPeak > sData.imdata(i).roiMeta(roi).trackAlignedPeak
        lickCoupledROIs(i) =  roi;
        i = i + 1;
    end
    if  sData.imdata(i).roiMeta(roi).activeTrialFractionInBlock > 0.3 && sData.imdata(i).roiMeta(roi).lickAlCtrlSMIDff > 3  && sData.imdata(i).roiMeta(roi).lickAlignedPeakDff > sData.imdata(i).roiMeta(roi).trackAlignedPeakDff
        lickCoupledROIsDff(j) =  roi;
        j = j + 1;
    end
    %}
    
    
end


% Analyze ROIs

% Active ROIs
threshold = 0.3;
%freqThreshold = 1.5;

% activeROIs = intersect(find(activeTrialFractionCtrl > threshold),find(peakFreqsCtrl > freqThreshold));
activeROIs = find(activeTrialFraction > threshold);
inactiveROIs = setdiff(1:nROIs,activeROIs)';
maxSMI = max(SMIDeconv');
tunedROIs = find(maxSMI > 3);
% tunedROIs = find(lickAlCtrlSMI > 3);
untunedROIs = setdiff(1:nROIs,tunedROIs)';



%% Plot pie charts
%{
data = NaN;
data(1) = numel(intersect(inactiveROIs,untunedROIs));
data(2) = numel(intersect(inactiveROIs,tunedROIs));
data(3) = numel(intersect(activeROIs,tunedROIs));
data(4) = numel(intersect(activeROIs,untunedROIs));

numel(intersect(find(activeTrialFractionCtrl > threshold),find(peakFreqsCtrl > freqThreshold)))

label = {'Inactive untuned','Inactive tuned','Active tuned','Active untuned'};


figure
subplot(2,2,1)
pie(data,[1 0 0 0],label)
subplot(2,2,3)
pie(data,[0 1 0 0],label)
subplot(2,2,4)
pie(data,[0 0 1 0],label)
subplot(2,2,2)
pie(data,[0 0 0 1],label)


figure
plot(activeTrialFractionCtrl,peakFreqs,'.')

%}
% numel(find(activeTrialFractionCtrl>0.3))
% numel(find(peakFreqs>1))
% numel(intersect(find(activeTrialFractionCtrl>0.3),find(peakFreqs'>1)))


%threshold = 0.3;
%activeROIs = find(activeTrialFraction > threshold);
%nactiveROIs = numel(activeROIs)/nROIs;

%figure; histogram(activeTrialFraction)

%% Save data to sData

% save(fullfile(sDataDir,sessionID),'sData');

% sData.imdata(i).roiMeta = roiMeta;

sData.imdata(i).binnedRoisDeconv = binnedRoisDeconv;
sData.imdata(i).binnedRoisDeconvRate = binnedRoisDeconvRate;
sData.imdata(i).binnedRoisDff = binnedRoisDff;
sData.imdata(i).binNumber = binNumber;
sData.imdata(i).binSize = binSize;
sData.imdata(i).nAllTrials = allTrials;

meta.fieldWidthThreshold = fieldWidthThreshold;
meta.firstThreshold = firstThreshold;
meta.reliabilityThreshold = reliabilityThreshold;
meta.smoothSpan = smoothSpan;
sData.imdata(i).placeCells.meta = meta;

%sData.imdata(i).placeCells.nPlaceCells = nPlaceCells;
%sData.imdata(i).placeCells.placeCellFraction = placeCellFraction;
%sData.imdata(i).placeCells.placeCells = placeCells;
% sData.imdata(i).placeCells.sortedPlaceCells = sortedPlaceCells;

sData.imdata(i).peakActPos = peakActPos;
sData.imdata(i).peakActPosInTrials = peakActPosInTrials;

% peakBinDeconv
% lickAlPeakBin
% peakBinDff
% lickAlPeakBinDff

% ctrlSMIDff
% lickAlCtrlSMIDff
sData.imdata(i).activeROIs = activeROIs;
sData.imdata(i).inactiveROIs = inactiveROIs;
sData.imdata(i).tunedROIs = tunedROIs;
sData.imdata(i).untunedROIs = untunedROIs;


end

%}
