function data = bb_MVPA_Analysis(study, data, subNum, activation)

%%%%%%%%%%%%%%%%%%%
% MVPA by BV
%%%%%%%%%%%%%%%%%%%
if isfield(data.sub(subNum).study.stats, 'bv') && isfield(study.stats, 'bv')
    for ii=1:length(data.sub(subNum).study.stats.bv)
        fprintf('Running %s bv(%d)\n', data.sub(subNum).info.name, ii);
        data.sub(subNum).study.stats.bv = mainfx(study, data, data.sub(subNum).study.stats.bv, ii, subNum, study.stats.bv(ii), activation);
    end
end

%%%%%%%%%%%%%%%%%%%
% MVPA by IV
%%%%%%%%%%%%%%%%%%%
if isfield(data.sub(subNum).study.stats, 'iv') && isfield(study.stats, 'iv')
    for ii=1:length(data.sub(subNum).study.stats.iv)
        fprintf('Running %s iv(%d)\n', data.sub(subNum).info.name, ii);
        data.sub(subNum).study.stats.iv = mainfx(study, data, data.sub(subNum).study.stats.iv, ii, subNum, study.stats.iv(ii), activation);
    end
end

%%%%%%%%%%%%%%%%%%%
% MVPA by DV
%%%%%%%%%%%%%%%%%%%
if isfield(data.sub(subNum).study.stats, 'dv') && isfield(study.stats, 'dv')
    for ii=1:length(data.sub(subNum).study.stats.dv)
        fprintf('Running %s dv(%d)\n', data.sub(subNum).info.name, ii);
        data.sub(subNum).study.stats.dv = mainfx(study, data, data.sub(subNum).study.stats.dv, ii, subNum, study.stats.dv(ii), activation);
    end
end

%%%%%%%%%%%%%%%%%%%
% MVPA by MV
%%%%%%%%%%%%%%%%%%%
if isfield(data.sub(subNum).study.stats, 'mv') && isfield(study.stats, 'mv')
    for ii=1:length(data.sub(subNum).study.stats.mv)
        fprintf('Running %s mv(%d)\n', data.sub(subNum).info.name, ii);
        data.sub(subNum).study.stats.mv = mainfx(study, data, data.sub(subNum).study.stats.mv, ii, subNum, study.stats.mv(ii), activation);
    end
end
end

%%
function datavar = mainfx(study, data, datavar, ii, subNum, studyvar, activation)

folder = sprintf('%s/%s/mvpa_results/%s', data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, studyvar.tag);
if exist(folder, 'dir')~=7
    mkdir(folder);
end

% Analyze Levels
fprintf('Analyzing levels');
if isfield(datavar(ii), 'level')
    for ilevel=1:length(datavar(ii).level)
        if isfield(datavar(ii).level(ilevel), 'types')
            clear temp;
            for tt = 1:length(datavar(ii).level(ilevel).types)
                temp(tt,:) = datavar(ii).level(ilevel).types(tt) == data.sub(subNum).cond_sequence;
            end
            diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
            if exist('temp') && sum(sum(temp)) >= 4
                blocks = ismember(data.sub(subNum).cond_sequence, datavar(ii).level(ilevel).types);
                matrix1 = activation.session.stim_act(blocks,:);
                matrix2 = activation.session.null_act(blocks,:);
                datavar(ii).level(ilevel).mvpa.temp = [];
                %fprintf('%s  Running bb_MVPA_ROI_Analysis on %s Level %d\n', datestr(now), data.sub(subNum).info.name, ilevel);
                %[datavar, datavar(ii).level(ilevel).mvpa] = bb_MVPA_ROI_Analysis(study, data, subNum, datavar, ii, datavar(ii).level(ilevel).mvpa, matrix1, matrix2, []);
                fprintf('%s  Running bb_MVPA_Searchlight_Analysis on %s Level %d\n', datestr(now), data.sub(subNum).info.name, ilevel);
                datavar(ii).level(ilevel).mvpa = bb_MVPA_Searchlight_Analysis(data, subNum, datavar(ii).tag, datavar(ii).level(ilevel).mvpa, 0.0001, datavar(ii).level(ilevel).name, matrix1, matrix2);
                clear matrix1; clear matrix2;
            end
        end
    end
end

% Analyze contrasts
fprintf('Analyzing contrasts');
if isfield(datavar(ii), 'contrast')
    % Generate useful data about each contrast
    fprintf('%s  Running bb_MVPA_Contrast_Prep on %s', datestr(now), data.sub(subNum).info.name);
    datavar = bb_MVPA_Contrast_Prep(data, subNum, datavar, ii);
    %Create a template for the confusion matrices
    fprintf('%s  Running bb_MVPA_Confusion_Matrix_Prep on %s', datestr(now), data.sub(subNum).info.name);
    datavar = bb_MVPA_Confusion_Matrix_Prep(study, data, subNum, datavar, ii);
    
    if isfield(datavar(ii), 'mvpa') && isfield(datavar(ii).mvpa, 'confusion_matrix')
        for icontrast=1:length(datavar(ii).contrast)
            if datavar(ii).contrast(icontrast).mvpa.suffTrials
                matrix1 = activation.session.betas(datavar(ii).contrast(icontrast).mvpa.cond1_blocks,:); %Create a matrix containing only condition 1 betas.
                matrix2 = activation.session.betas(datavar(ii).contrast(icontrast).mvpa.cond2_blocks,:); %Create a matrix containing only condition 2 betas.
                fprintf('%s  Running bb_MVPA_ROI_Analysis on %s Contrast %d', datestr(now), data.sub(subNum).info.name, icontrast);
                [datavar, datavar(ii).contrast(icontrast).mvpa] = bb_MVPA_ROI_Analysis(study, data, subNum, datavar, ii, datavar(ii).contrast(icontrast).mvpa, matrix1, matrix2, datavar(ii).contrast(icontrast).mvpa.conf_mat_loc);
                fprintf('%s  Running bb_MVPA_Searchlight_Analysis on %s Contrast %d', datestr(now), data.sub(subNum).info.name, icontrast);
                datavar(ii).contrast(icontrast).mvpa = bb_MVPA_Searchlight_Analysis(data, subNum, datavar(ii).tag, datavar(ii).contrast(icontrast).mvpa, 0.0001, sprintf('Contrast_%d', icontrast), matrix1, matrix2);
                clear matrix1; clear matrix2;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%REMOVING FOR TESTING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %Fill in spots on confusion matrix that are not an official "contrast"
%         for icond1=1:length(datavar(ii).mvpa.confusion_matrix.cond)
%             for icond2=1:(length(datavar(ii).mvpa.confusion_matrix.cond)-icond1)
%                 if isempty(datavar(ii).mvpa.confusion_matrix.contrast_tracker(icond1, icond2).contrast)
%                     matrix1 = activation.session.betas(datavar(ii).mvpa.confusion_matrix.cond(icond1).blocks,:);
%                     matrix2 = activation.session.betas(datavar(ii).mvpa.confusion_matrix.cond(icond2).blocks,:);
%                     [datavar, ~] = bb_MVPA_ROI_Analysis(data, subNum, datavar, ii, [], matrix1, matrix2, [icond1,icond2]);
%                     clear matrix1; clear matrix2;
%                 end
%             end
%         end
        %Generate a confusion matrix figure
        fprintf('%s  Running bb_MVPA_Conf_Mat_Figure on %s', datestr(now), data.sub(subNum).info.name);
        datavar = bb_MVPA_Conf_Mat_Figure(study, data, subNum, datavar, ii);
    end
end
end
