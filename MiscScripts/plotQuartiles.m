
function varargout = plotQuartiles(x,trialMatrix,color,smoothingWindow,lineWidth)
% based on: http://jvoigts.scripts.mit.edu/blog/nice-shaded-plots/
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

quartiles = quantile(trialMatrix,[0.25 0.50 0.75],1);
quartiles = smoothdata(quartiles,2,'gaussian',smoothingWindow,'omitnan');
med = quartiles(2,:);

patchX = [x, fliplr(x)];
patchY = [quartiles(1,:), fliplr(quartiles(3,:))];

patchX = patchX(~isnan(patchY));
patchY = patchY(~isnan(patchY));

p = patch(patchX,patchY,1,'FaceColor',color,'EdgeColor','none');
p.FaceAlpha = 0.2; % Set transparency

else
med = smoothdata(trialMatrix,'gaussian',smoothingWindow,'omitnan');

end

hold on
plot(x,med,'Color',color,'LineWidth',lineWidth);


    


end