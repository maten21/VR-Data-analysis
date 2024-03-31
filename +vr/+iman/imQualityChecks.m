function sData = imQualityChecks(sData, sDataDir)

%% quality check

i = 1; 


nFOVs = length(sData.imdata);
originalSessionID = sData.sessionInfo.sessionID;

sData2 = sData;
if nFOVs > 1
    for i = 1:1:nFOVs
        sData2.imdata = sData.imdata(i);
        sData2.sessionInfo.sessionID = [originalSessionID(1:17), '-FOV', num2str(i)];
        if ~isfield(sData2.imdata, 'roiStat')
            roiStat = getRoiActivityStats(sData2);
            sData2.imdata.roiStat = roiStat;
            sData.imdata(i).roiStat = roiStat;
        end

        
        plotRoiActivityStats(sData2)
        
        if nargin == 2
            saveas(gcf,fullfile(sDataDir,['roiStats', '_fov', num2str(i), '.png']));
        end
    end
    
else

        sData2.sessionInfo.sessionID = originalSessionID(1:17);
        if ~isfield(sData2.imdata, 'roiStat')
            roiStat = getRoiActivityStats(sData2);
            sData2.imdata.roiStat = roiStat;
            sData.imdata.roiStat = roiStat;              
        end
        
        
        plotRoiActivityStats(sData2)
        
        if nargin == 2
            saveas(gcf,fullfile(sDataDir,['roiStats', '_fov', num2str(i), '.png']));
        end


end

end

function plotRoiActivityStats(sData)
%plotRoiActivityStats Make figure with roi statistics for a session
%
%   plotRoiActivityStats(sData) plots distributions of the different
%   parameters in violin plots.
%
%   See also getRoiActivityStats


    f = openfig('roiActivityStatsSummary.fig');
    hAx = findobj(f, 'Type', 'Axes');
    
    if isfield(sData.imdata, 'roiStat')
        S = sData.imdata.roiStat;
    else
        S = getRoiActivityStats(sData);
    end

    sTitle = sData.sessionInfo.sessionID;
    figtitle(f, sTitle, 16, 0.97)

    hVio = gobjects(3,1);
    [hVio(1)] = violin(hAx(1), S.peakDff, 'mc', [], 'plotlegend', 0);
    [hVio(2)] = violin(hAx(2), S.signalToNoise, 'mc', [], 'plotlegend', 0);
    [hVio(3)] = violin(hAx(3), S.activityLevel, 'mc', [], 'plotlegend', 1);
    
    cmap = viridis();
    for i = 1:3
        c = cmap(randi([1,256], 1), :);
        set(hVio(i), 'FaceColor', c, 'EdgeColor', c)
    end
    
    set(hAx, 'XTick', [], 'XLim', [0.6, 1.4])
    
    hAx(1).Title.String = 'Peak DFF';
    hAx(2).Title.String = 'SNR';
    hAx(3).Title.String = 'tActive/tTotal';
    
    set([hAx.Title], 'FontWeight', 'normal')


end