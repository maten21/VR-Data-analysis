% This finction filters the raw lick signal and creates a summary figure of the lick signal data
% Input: raw lickSignal, samplingRate, lengthThreshold, freqThreshold   
% Output: lickEvent signal, lickCount, licErrorPercentage, lickFig


function [lickEvents,lickCount,licErrorPercentage,lickFig] = analyseLickData(lickSignal,samplingRate,lickLengthThreshold,lickFreqThreshold)

lickStartIndexesRaw = find(diff(lickSignal) == 1)+1;
lickEndIndexesRaw = find(diff(lickSignal) == -1);

% Correct for licks not fully captured at the beginning or end of the recording
if lickSignal(1) == 1
    lickEndIndexesRaw = lickEndIndexesRaw(2:numel(lickEndIndexesRaw));
end
if lickSignal(numel(lickSignal)) == 1
    lickStartIndexesRaw = lickStartIndexesRaw(1:numel(lickStartIndexesRaw)-1);
end

lickLengthsMsRaw = (lickEndIndexesRaw - lickStartIndexesRaw)/(samplingRate/1000); %Convert data to ms

diffLickStart = diff(lickStartIndexesRaw);
lickFreqRaw = samplingRate./diffLickStart;

shortLickIndexes = lickStartIndexesRaw(find(lickLengthsMsRaw < lickLengthThreshold));
tooFastLickIndexes = lickStartIndexesRaw(find(lickFreqRaw > lickFreqThreshold));
lickErrorsIndexes = union(shortLickIndexes,tooFastLickIndexes);

% Filtered lick start indexes
lickStartIndexes = setdiff(lickStartIndexesRaw,lickErrorsIndexes);

lickCount = numel(lickStartIndexes);
lickErrors = numel(lickErrorsIndexes);
licErrorPercentage = lickErrors/numel(lickStartIndexesRaw)*100;
% Create lick event array
lickEvents(1:numel(lickSignal)) = zeros;
lickEvents(lickStartIndexes) = 1;


try
% PLOT LICK ANALYSIS SUMMARY FIGURE
% Plot lick length distribution
lickFig = figure('Color','white','Position',[0 0 1200 750]);
subplot(2,3,1);
histogram(lickLengthsMsRaw)
xlim([0 100]);
xlabel('Lick length (ms)');
ylabel('Count in bin');

% Check wether lick lengths determine with what frequency the next lick follows
subplot(2,3,4);
rectangle('Position',[lickLengthThreshold,0.5,0.5,1000],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
hold on
rectangle('Position',[0.5,lickFreqThreshold,100,0.5],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
hold on
plot(lickLengthsMsRaw(1:numel(lickFreqRaw)),lickFreqRaw,'.')
dim = [0.24 0.14 0.3 0.3];
str = {['Lick count: ', num2str(lickCount)],['Lick errors: ', num2str(lickErrors)],[num2str(licErrorPercentage), ' % errors']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

xlim([0 100]);
ylim([0 100]);
xlabel('Lick length (ms)');
ylabel('Instantaneous frequency (Hz)');


% Plot lickSignal aligned to short licks
shortLicksMatrix = NaN(lickErrors,1500);
for i = 1:1:lickErrors
    shortLicksMatrix(i,:) = lickSignal(lickErrorsIndexes(i)-599:lickErrorsIndexes(i)+900);
end

subplot(2,3,[2 3 5 6]);
mymap = [1 1 1 ; lines(1)];
imagesc(-199:1:300,1:1:numel(lickLengthsMsRaw(lickLengthsMsRaw<lickLengthThreshold)),shortLicksMatrix)
colormap(gca,mymap);
colormap(gca,mymap);
xlabel('Lick length (ms)');
ylabel('Erroneous licks');

end

%{

% PLOT LICK ANALYSIS SUMMARY FIGURE
% Plot lick length distribution
figure('Color','white','Position',[0 0 1200 750]);
subplot(2,3,4);
histogram(lickLengthsMs)
xlim([0 100]);

% Check wether lick lengths determine with what frequency the next lick follows
subplot(2,3,4);
rectangle('Position',[lickLengthThreshold,0.5,0.5,1000],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
hold on
rectangle('Position',[0.5,lickFreqThreshold,100,0.5],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone
hold on
plot(lickLengthsMsRaw(1:numel(lickFreqRaw)),lickFreqRaw,'.')
xlim([0 100]);
ylim([0 100]);

% Plot lickSignal aligned to short licks
shortLicksMatrix = NaN(lickErrors,1500);
for i = 1:1:lickErrors
    shortLicksMatrix(i,:) = lickSignal(shortLickIndexes(i)-599:shortLickIndexes(i)+900);
end

subplot(2,3,[2 3 5 6]);
mymap = [1 1 1 ; lines(1)];
imagesc(-199:1:300,1:1:numel(lickLengthsMsRaw(lickLengthsMsRaw<lickLengthThreshold)),shortLicksMatrix)
colormap(gca,mymap);
%}
end