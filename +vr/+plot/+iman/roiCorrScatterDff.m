function [] = roiCorrScatterDff(sDataFiles,filePath)

if nargin < 2
    clear
    [sDataFiles, filePath] = vr.loadData('light');
end

if ~iscell(sDataFiles)
    sDataFiles = {sDataFiles};
end


for f = 1:1:length(sDataFiles)
    
    for fov = 1:1:length(sDataFiles{1, f}.imdata)
       
        %A = nan(sDataFiles{1, f}.imdata(fov).nROIs,1);
        %AB = nan(sDataFiles{1, f}.imdata(fov).nROIs,1);
        %B = nan(sDataFiles{1, f}.imdata(fov).nROIs,1);
        
        %for roi = 1:1:sDataFiles{1, f}.imdata(fov).nROIs
            A =      [sDataFiles{1, f}.imdata(fov).roiMeta.identPartCorrCoefA];
            B =      [sDataFiles{1, f}.imdata(fov).roiMeta.identPartCorrCoefB];
            AB =     [sDataFiles{1, f}.imdata(fov).roiMeta.identPartCorrCoefAB];
            ASign =  [sDataFiles{1, f}.imdata(fov).roiMeta.identPartIsSignCorrA];
            BSign =  [sDataFiles{1, f}.imdata(fov).roiMeta.identPartIsSignCorrB];
            ABSign = [sDataFiles{1, f}.imdata(fov).roiMeta.identPartIsSignCorrAB];
        %end
        
        
            ASignRois = find(ASign);
            BSignRois = find(BSign);
            ABSignRois = find(ABSign);
            nROIs = sDataFiles{1, f}.imdata(fov).nROIs;
            
            remapRoisA = setdiff(ASignRois,ABSignRois); 
            nonRemapRoisA = intersect(ASignRois,ABSignRois);
            remapRoisB = setdiff(BSignRois,ABSignRois); 
            nonRemapRoisB = intersect(BSignRois,ABSignRois);
            
            tunedFractionA = numel(ASignRois) / nROIs;
            untunedFractionA = 1 - tunedFractionA;
            remappingFractionA = numel(setdiff(ASignRois,ABSignRois)) / numel(ASignRois);
            tunedFractionB = numel(BSignRois) / nROIs;
            untunedFractionB = 1 - tunedFractionB;
            remappingFractionB = numel(setdiff(BSignRois,ABSignRois)) / numel(BSignRois);
            
            
            
            
            % only significant in AB probably these are low activity ROIs
            % so even odd does not correlate

        
        barData = NaN;
        
        barData(1,1) = untunedFractionA;
        barData(1,2) = numel(intersect(ASignRois,ABSignRois)) / nROIs; % non remapping tuned rois
        barData(1,3) = numel(setdiff(ASignRois,ABSignRois)) / nROIs; % remapping tuned rois

        barData(2,1) = untunedFractionB;
        barData(2,2) = numel(intersect(BSignRois,ABSignRois)) / nROIs; % non remapping tuned rois
        barData(2,3) = numel(setdiff(BSignRois,ABSignRois)) / nROIs; % remapping tuned rois
        
        
       
        
        
        figure('Color','white','Position',[0 0 900 300])
                               
        h = subplot(1,3,1);
        hold on
        h.PlotBoxAspectRatio = [1 1 1];
        plot(A,B,'.K')
        plot(A(ABSignRois),B(ABSignRois),'.')
        plot([A(remapRoisA), A(remapRoisB)],[B(remapRoisA), B(remapRoisB)],'.')
        plot([-1 1],[-1 1],'--k')
        xlim([-1 1])
        ylim([-1 1])
        xlabel('Corr. A (odd x even)')
        ylabel('Corr. B (odd x even)')
        title('Same vs. same')
        
       
        h = subplot(1,3,2);
        hold on
        h.PlotBoxAspectRatio = [1 1 1];
        plot([A(nonRemapRoisA), B(nonRemapRoisB)],[AB(nonRemapRoisA), AB(nonRemapRoisB)],'.')
        plot([A(remapRoisA), B(remapRoisB)],[AB(remapRoisA), AB(remapRoisB)],'.')
%        plot(A(ABSign),AB(ABSign),'.')
        plot([0 1],[0 1],'--k')
        plot([0 1],[0 0],'-k')
        xlim([0.5 1])
        %ylim([-1 1])
        xlabel('Corr. coef. AxA or BxB')
        ylabel('Corr. coef. AxB')
        title('Same vs. different')
        
        
        h = subplot(1,3,3);
        hold on
        h.PlotBoxAspectRatio = [1 1 1];
        
        x = categorical({'Familiar (A)','New (B)'});
        bar(x,barData(:,2:3),'stacked')
        ylim([0 1])
        ylabel('Fraction of ROIs')
        xlabel('Tuned ROIs Trial blocks')
        legend({'Stable (AxB corr.)','Remap (AxB not corr.)'},'Location','northeast')
        title('Fraction of neurons')

        suptitle(['\fontsize{12}' sDataFiles{1, f}.imdata(fov).fovLocation ' - ' sDataFiles{1, f}.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' - Identical part' ])
     

        
        if ~isfolder(fullfile(filePath,'RoiCorrScatterDff')); mkdir(fullfile(filePath,'RoiCorrScatterDff')); end
        saveas(gcf,strcat(fullfile(filePath,'RoiCorrScatterDff',['fov' num2str(fov) '-' sDataFiles{1, f}.imdata(fov).fovLocation '-' 'RoiCorrScatterIdentPart_' sDataFiles{1, f}.sessionInfo.sessionID(1:17)]),'.png'));
        
        close gcf
        
        
        
        
        
        
        
        
         %A = nan(sDataFiles{1, f}.imdata(fov).nROIs,1);
        %AB = nan(sDataFiles{1, f}.imdata(fov).nROIs,1);
        %B = nan(sDataFiles{1, f}.imdata(fov).nROIs,1);
        
        %for roi = 1:1:sDataFiles{1, f}.imdata(fov).nROIs
            A =      [sDataFiles{1, f}.imdata(fov).roiMeta.uniquePartCorrCoefA];
            B =      [sDataFiles{1, f}.imdata(fov).roiMeta.uniquePartCorrCoefB];
            AB =     [sDataFiles{1, f}.imdata(fov).roiMeta.uniquePartCorrCoefAB];
            ASign =  [sDataFiles{1, f}.imdata(fov).roiMeta.uniquePartIsSignCorrA];
            BSign =  [sDataFiles{1, f}.imdata(fov).roiMeta.uniquePartIsSignCorrB];
            ABSign = [sDataFiles{1, f}.imdata(fov).roiMeta.uniquePartIsSignCorrAB];
        %end
        
        
            ASignRois = find(ASign);
            BSignRois = find(BSign);
            ABSignRois = find(ABSign);
            nROIs = sDataFiles{1, f}.imdata(fov).nROIs;
            
            remapRoisA = setdiff(ASignRois,ABSignRois); 
            nonRemapRoisA = intersect(ASignRois,ABSignRois);
            remapRoisB = setdiff(BSignRois,ABSignRois); 
            nonRemapRoisB = intersect(BSignRois,ABSignRois);
            
            tunedFractionA = numel(ASignRois) / nROIs;
            untunedFractionA = 1 - tunedFractionA;
            remappingFractionA = numel(setdiff(ASignRois,ABSignRois)) / numel(ASignRois);
            tunedFractionB = numel(BSignRois) / nROIs;
            untunedFractionB = 1 - tunedFractionB;
            remappingFractionB = numel(setdiff(BSignRois,ABSignRois)) / numel(BSignRois);
            
            
            
            
            % only significant in AB probably these are low activity ROIs
            % so even odd does not correlate

        
        barData = NaN;
        
        barData(1,1) = untunedFractionA;
        barData(1,2) = numel(intersect(ASignRois,ABSignRois)) / nROIs; % non remapping tuned rois
        barData(1,3) = numel(setdiff(ASignRois,ABSignRois)) / nROIs; % remapping tuned rois

        barData(2,1) = untunedFractionB;
        barData(2,2) = numel(intersect(BSignRois,ABSignRois)) / nROIs; % non remapping tuned rois
        barData(2,3) = numel(setdiff(BSignRois,ABSignRois)) / nROIs; % remapping tuned rois
        
        
       
        
        
        figure('Color','white','Position',[0 0 900 300])
                               
        h = subplot(1,3,1);
        hold on
        h.PlotBoxAspectRatio = [1 1 1];
        plot(A,B,'.K')
        plot(A(ABSignRois),B(ABSignRois),'.')
        plot([A(remapRoisA), A(remapRoisB)],[B(remapRoisA), B(remapRoisB)],'.')
        plot([-1 1],[-1 1],'--k')
        xlim([-1 1])
        ylim([-1 1])
        xlabel('Corr. A (odd x even)')
        ylabel('Corr. B (odd x even)')
        title('Same vs. same')
        
       
        h = subplot(1,3,2);
        hold on
        h.PlotBoxAspectRatio = [1 1 1];
        plot([A(nonRemapRoisA), B(nonRemapRoisB)],[AB(nonRemapRoisA), AB(nonRemapRoisB)],'.')
        plot([A(remapRoisA), B(remapRoisB)],[AB(remapRoisA), AB(remapRoisB)],'.')
%        plot(A(ABSign),AB(ABSign),'.')
        plot([0 1],[0 1],'--k')
        plot([0 1],[0 0],'-k')
        xlim([0.5 1])
        %ylim([-1 1])
        xlabel('Corr. coef. AxA or BxB')
        ylabel('Corr. coef. AxB')
        title('Same vs. different')
        
        
        h = subplot(1,3,3);
        hold on
        h.PlotBoxAspectRatio = [1 1 1];
        
        x = categorical({'Familiar (A)','New (B)'});
        bar(x,barData(:,2:3),'stacked')
        ylim([0 1])
        ylabel('Fraction of ROIs')
        xlabel('Tuned ROIs Trial blocks')
        legend({'Stable (AxB corr.)','Remap (AxB not corr.)'},'Location','northeast')
        title('Fraction of neurons')

        suptitle(['\fontsize{12}' sDataFiles{1, f}.imdata(fov).fovLocation ' - ' sDataFiles{1, f}.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' - Unique part' ])
     

        
        if ~isfolder(fullfile(filePath,'RoiCorrScatterDff')); mkdir(fullfile(filePath,'RoiCorrScatterDff')); end
        saveas(gcf,strcat(fullfile(filePath,'RoiCorrScatterDff',['fov' num2str(fov) '-' sDataFiles{1, f}.imdata(fov).fovLocation '-' 'RoiCorrScatterUniquePart_' sDataFiles{1, f}.sessionInfo.sessionID(1:17)]),'.png'));
        
        close gcf
        
        
    end
end

end