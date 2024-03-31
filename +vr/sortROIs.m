function [sortedRois, sortIndexes] = sortROIs(posTunCurves,subsetOfInterest,smoothSpan)
%% Sort ROIs based on peak activity position of the tuning curves
%
% Required inputs: 
%   - posTunCurves matrix (nROIs X nSpatialBins) 
% 
% Optional inputs: 
%   - smoothSpan (default = 1)
%   - subsetOfInterest (e.g. place cells or active ROIs)
%
% updated 2022.06.20

%if ~exist('posTunCurves') % default deconv light off trials  
%    posTunCurves = sData.imdata.avBinnedRois.avBinnedRoisBCNLightOffDeconv;
%end

if exist('smoothSpan') % Initial smoothing of tuning curves
    posTunCurves =  smoothdata(posTunCurves,2,'gaussian',smoothSpan);
end
%if exist('subsetOfInterest') % apply subset
%    posTunCurves = posTunCurves(subsetOfInterest,:);    
%end



nROIs = size(posTunCurves,1);

% determine place field peaks:
[~, peakActPos] = max(posTunCurves'); % [maxVal, index] = max(A)


sortedRois(nROIs,2) = zeros;
sortedRois(:,2) = peakActPos;
sortedRois(:,1) = 1:1:nROIs;

[sortedRois, sortIndexes] = sortrows(sortedRois,2);
sortedRois = sortedRois(:,1);


if exist('subsetOfInterest') % apply subset
    sortedRois = sortedRois(find(ismember(sortedRois,subsetOfInterest)));
    sortIndexes = sortIndexes(find(ismember(sortedRois,subsetOfInterest)));
end


end