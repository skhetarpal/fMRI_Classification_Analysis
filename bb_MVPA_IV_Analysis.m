function data = bb_MVPA_IV_Analysis(study, data, subNum, activation)

if isfield(data.sub(subNum).study.stats, 'iv') && isfield(study.stats, 'iv')
    for ii=1:length(data.sub(subNum).study.stats.iv)
        fprintf('Running %s iv(%d)\n', data.sub(subNum).info.name, ii);
        %data.sub(subNum).study.stats.iv = mainfx(study, data, data.sub(subNum).study.stats.iv, ii, subNum, study.stats.iv(ii), activation);
        %function data.sub(subNum).study.stats.iv = mainfx(study, data, datavar, ii, subNum, studyvar, activation)
        
        folder = sprintf('%s/%s/mvpa_results/%s', data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, study.stats.iv(ii).tag);
        if exist(folder, 'dir')~=7
            mkdir(folder);
        end
        
        % Analyze Levels
        fprintf('Analyzing levels');
        if isfield(data.sub(subNum).study.stats.iv(ii), 'level')
            for ilevel=1:length(data.sub(subNum).study.stats.iv(ii).level)
                if isfield(data.sub(subNum).study.stats.iv(ii).level(ilevel), 'types')
                    clear temp;
                    for tt = 1:length(data.sub(subNum).study.stats.iv(ii).level(ilevel).types)
                        temp(tt,:) = data.sub(subNum).study.stats.iv(ii).level(ilevel).types(tt) == data.sub(subNum).cond_sequence;
                    end
                    diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
                    if exist('temp') && sum(sum(temp)) >= 4
                        blocks = ismember(data.sub(subNum).cond_sequence, data.sub(subNum).study.stats.iv(ii).level(ilevel).types);
                        matrix1 = activation.session.stim_act(blocks,:);
                        matrix2 = activation.session.null_act(blocks,:);
                        data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa.temp = [];
                        if ~isfield(data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa, 'roi_done') || data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa.roi_done ~= 0
                            fprintf('%s  Running bb_MVPA_ROI_Analysis on %s Level %d\n', datestr(now), data.sub(subNum).info.name, ilevel);
                            [data.sub(subNum).study.stats.iv, data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa] = ...
                                bb_MVPA_ROI_Analysis(study, data, subNum, data.sub(subNum).study.stats.iv, ii, ...
                                data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa, matrix1, matrix2, []);
                            data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa.roi_done = 1;
                            datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
                            save(datafile, 'data');
                        end
                        if ~isfield(data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa, 'sl_done') || data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa.sl_done ~= 0
                            fprintf('%s  Running bb_MVPA_Searchlight_Analysis on %s Level %d\n', datestr(now), data.sub(subNum).info.name, ilevel);
                            data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa = ...
                                bb_MVPA_Searchlight_Analysis(data, subNum, data.sub(subNum).study.stats.iv(ii).tag, ...
                                data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa, 0.0001, ...
                                data.sub(subNum).study.stats.iv(ii).level(ilevel).name, matrix1, matrix2);
                            data.sub(subNum).study.stats.iv(ii).level(ilevel).mvpa.sl_done = 1;
                            datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
                            save(datafile, 'data');
                        end
                        clear matrix1; clear matrix2;
                    end
                end
            end
        end
        
        % Analyze contrasts
        fprintf('Analyzing contrasts');
        if isfield(data.sub(subNum).study.stats.iv(ii), 'contrast')
            % Generate useful data about each contrast
            fprintf('%s  Running bb_MVPA_Contrast_Prep on %s', datestr(now), data.sub(subNum).info.name);
            data.sub(subNum).study.stats.iv = bb_MVPA_Contrast_Prep(data, subNum, data.sub(subNum).study.stats.iv, ii);
            %Create a template for the confusion matrices
            fprintf('%s  Running bb_MVPA_Confusion_Matrix_Prep on %s', datestr(now), data.sub(subNum).info.name);
            data.sub(subNum).study.stats.iv = bb_MVPA_Confusion_Matrix_Prep(study, data, subNum, data.sub(subNum).study.stats.iv, ii);
            
            if isfield(data.sub(subNum).study.stats.iv(ii), 'mvpa') && isfield(data.sub(subNum).study.stats.iv(ii).mvpa, 'confusion_matrix')
                for icontrast=1:length(data.sub(subNum).study.stats.iv(ii).contrast)
                    if data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.suffTrials
                        matrix1 = activation.session.betas(data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.cond1_blocks,:); %Create a matrix containing only condition 1 betas.
                        matrix2 = activation.session.betas(data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.cond2_blocks,:); %Create a matrix containing only condition 2 betas.
                        if ~isfield(data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa, 'roi_done') || data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.roi_done ~= 0
                            fprintf('%s  Running bb_MVPA_ROI_Analysis on %s Contrast %d', datestr(now), data.sub(subNum).info.name, icontrast);
                            [data.sub(subNum).study.stats.iv, data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa] = ...
                                bb_MVPA_ROI_Analysis(study, data, subNum, data.sub(subNum).study.stats.iv, ii, ...
                                data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa, matrix1, matrix2, ...
                                data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.conf_mat_loc);
                            data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.roi_done = 1;
                            datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
                            save(datafile, 'data');
                        end
                        if ~isfield(data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa, 'sl_done') || data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.sl_done ~= 0
                            fprintf('%s  Running bb_MVPA_Searchlight_Analysis on %s Contrast %d', datestr(now), data.sub(subNum).info.name, icontrast);
                            data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa = ...
                                bb_MVPA_Searchlight_Analysis(data, subNum, data.sub(subNum).study.stats.iv(ii).tag, ...
                                data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa, 0.0001, sprintf('Contrast_%d', icontrast), ...
                                matrix1, matrix2);
                            data.sub(subNum).study.stats.iv(ii).contrast(icontrast).mvpa.sl_done = 1;
                            datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
                            save(datafile, 'data');
                        end
                        clear matrix1; clear matrix2;
                    end
                end
                
                %Fill in spots on confusion matrix that are not an official "contrast"
                for icond1=1:length(data.sub(subNum).study.stats.iv(ii).mvpa.confusion_matrix.cond)
                    for icond2=1:(length(data.sub(subNum).study.stats.iv(ii).mvpa.confusion_matrix.cond)-icond1)
                        if isempty(data.sub(subNum).study.stats.iv(ii).mvpa.confusion_matrix.contrast_tracker(icond1, icond2).contrast)
                            matrix1 = activation.session.betas(data.sub(subNum).study.stats.iv(ii).mvpa.confusion_matrix.cond(icond1).blocks,:);
                            matrix2 = activation.session.betas(data.sub(subNum).study.stats.iv(ii).mvpa.confusion_matrix.cond(icond2).blocks,:);
                            [data.sub(subNum).study.stats.iv, ~] = bb_MVPA_ROI_Analysis(data, subNum, data.sub(subNum).study.stats.iv, ii, [], matrix1, matrix2, [icond1,icond2]);
                            clear matrix1; clear matrix2;
                            datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
                            save(datafile, 'data');
                        end
                    end
                end
                %Generate a confusion matrix figure
                fprintf('%s  Running bb_MVPA_Conf_Mat_Figure on %s', datestr(now), data.sub(subNum).info.name);
                data.sub(subNum).study.stats.iv = bb_MVPA_Conf_Mat_Figure(study, data, subNum, data.sub(subNum).study.stats.iv, ii);
            end
        end
    end
end
end
