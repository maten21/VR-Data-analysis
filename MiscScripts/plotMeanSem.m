
function [pMean, pStd] = plotMeanSem(x,trialMatrix,color,smoothingWindow,lineWidth)
% http://jvoigts.scripts.mit.edu/blog/nice-shaded-plots/
% x: x coordinates
% y: either just one y vector, or 2xN or 3xN matrix of y-data
% fstr: format ('r' or 'b--' etc)
%
% example
% x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r');
 
if numel(x) ~= size(trialMatrix,2)
    trialMatrix = trialMatrix';
end
%{ 
switch nargin
    case 
    case
end
%}

if nargin < 2
    trialMatrix = x;
    x = 1:1:size(trialMatrix,2);
    color = lines(1);
    smoothingWindow = 1;
    lineWidth = 2;
end
if nargin < 3
    color = lines(1);
    smoothingWindow = 1;
    lineWidth = 2;
end
if nargin < 4
    smoothingWindow = 1;
    lineWidth = 2;
end
if nargin < 5
    lineWidth = 2;
end


if size(trialMatrix,2) > 1
    
    sessionMean = smoothdata(nanmean(trialMatrix),'gaussian',smoothingWindow,'omitnan');
    sessionStd = smoothdata(nanstd(trialMatrix),'gaussian',smoothingWindow,'omitnan');
    sessionSem = sessionStd / sqrt(height(trialMatrix));

    patchX = [x, fliplr(x)];
    patchY = [sessionMean + sessionSem, fliplr(sessionMean - sessionSem)];
    
    patchX = patchX(~isnan(patchY));
    patchY = patchY(~isnan(patchY));
    
    pStd = patch(patchX,patchY,1,'FaceColor',color,'EdgeColor','none');
    pStd.FaceAlpha = 0.2; % Set transparency
    
else
    sessionMean = smoothdata(trialMatrix,'gaussian',smoothingWindow,'omitnan');
end

hold on
pMean = plot(x,sessionMean,'Color',color,'LineWidth',lineWidth);



end