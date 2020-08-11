function datavar = bb_MVPA_Confusion_Matrix_Prep(study, data, subNum, datavar, ii)

% Build a master confusion matrix template from the contrasts
for icontrast=1:size(datavar(ii).contrast,2)
    if datavar(ii).contrast(icontrast).mvpa.suffTrials
        if ~isfield(datavar(ii), 'mvpa') || ~isfield(datavar(ii).mvpa, 'confusion_matrix'); %If this is the first contrast, set it to be the first and second conditions in the confusion matrix list.
            datavar(ii).mvpa.confusion_matrix.cond(1).name = datavar(ii).contrast(icontrast).mvpa.cond1_name;
            datavar(ii).mvpa.confusion_matrix.cond(2).name = datavar(ii).contrast(icontrast).mvpa.cond2_name;
            datavar(ii).mvpa.confusion_matrix.cond(1).levels = datavar(ii).contrast(icontrast).mvpa.cond1_levels;
            datavar(ii).mvpa.confusion_matrix.cond(2).levels = datavar(ii).contrast(icontrast).mvpa.cond2_levels;
            datavar(ii).mvpa.confusion_matrix.cond(1).blocks = datavar(ii).contrast(icontrast).mvpa.cond1_blocks;
            datavar(ii).mvpa.confusion_matrix.cond(2).blocks = datavar(ii).contrast(icontrast).mvpa.cond2_blocks;
            list_spot_1 = 1;
            list_spot_2 = 2;
        else %If this is not the first contrast, determine if the conditions are already on the list.  If not, add them.
            new_member1 = 1;
            new_member2 = 1;
            for icond = 1:length(datavar(ii).mvpa.confusion_matrix.cond)
                if strcmp(datavar(ii).contrast(icontrast).mvpa.cond1_name, datavar(ii).mvpa.confusion_matrix.cond(icond).name)
                    new_member1 = 0;
                    list_spot_1 = icond;
                end
                if strcmp(datavar(ii).contrast(icontrast).mvpa.cond2_name, datavar(ii).mvpa.confusion_matrix.cond(icond).name)
                    new_member2 = 0;
                    list_spot_2 = icond;
                end
            end
            if new_member1
                list_spot_1 = length(datavar(ii).mvpa.confusion_matrix.cond)+1;
                datavar(ii).mvpa.confusion_matrix.cond(list_spot_1).name = datavar(ii).contrast(icontrast).mvpa.cond1_name;
                datavar(ii).mvpa.confusion_matrix.cond(list_spot_1).levels = datavar(ii).contrast(icontrast).mvpa.cond1_levels;
                datavar(ii).mvpa.confusion_matrix.cond(list_spot_1).blocks = datavar(ii).contrast(icontrast).mvpa.cond1_blocks;
            end
            if new_member2
                list_spot_2 = length(datavar(ii).mvpa.confusion_matrix.cond)+1;
                datavar(ii).mvpa.confusion_matrix.cond(list_spot_2).name = datavar(ii).contrast(icontrast).mvpa.cond2_name;
                datavar(ii).mvpa.confusion_matrix.cond(list_spot_2).levels = datavar(ii).contrast(icontrast).mvpa.cond2_levels;
                datavar(ii).mvpa.confusion_matrix.cond(list_spot_2).blocks = datavar(ii).contrast(icontrast).mvpa.cond2_blocks;
            end
        end
        datavar(ii).mvpa.confusion_matrix.contrast_tracker(list_spot_1, list_spot_2).contrast = icontrast;
        datavar(ii).mvpa.confusion_matrix.contrast_tracker(list_spot_1, list_spot_2).order = 1;
        datavar(ii).mvpa.confusion_matrix.contrast_tracker(list_spot_2, list_spot_1).contrast = icontrast;
        datavar(ii).mvpa.confusion_matrix.contrast_tracker(list_spot_2, list_spot_1).order = 2;
        datavar(ii).contrast(icontrast).mvpa.conf_mat_loc = [list_spot_1, list_spot_2];
    end
end

% Initialize a confusion matrix for each roi.
if isfield(datavar(ii), 'mvpa') && isfield(datavar(ii).mvpa, 'confusion_matrix')
    for ap=1:length(study.part(data.sub(subNum).info.part).newapAnalysis)
        isess = study.part(data.sub(subNum).info.part).newapAnalysis(ap);
        for kk=1:length(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI)
            datavar(ii).mvpa.apAnalyses(isess).apROI{kk}.confusion_matrix.acc = NaN(length(datavar(ii).mvpa.confusion_matrix.cond),length(datavar(ii).mvpa.confusion_matrix.cond));
            datavar(ii).mvpa.apAnalyses(isess).apROI{kk}.confusion_matrix.pval = NaN(length(datavar(ii).mvpa.confusion_matrix.cond),length(datavar(ii).mvpa.confusion_matrix.cond));
        end
    end
end