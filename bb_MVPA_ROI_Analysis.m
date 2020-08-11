function [datavar, datavar2] = bb_MVPA_ROI_Analysis(study, data, subNum, datavar, ii, datavar2, matrix1, matrix2, conf_mat_loc)

for ap=1:length(study.part(data.sub(subNum).info.part).newapAnalysis)
    isess = study.part(data.sub(subNum).info.part).newapAnalysis(ap);
    for kk=1:length(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI)
        [foldera, subroota, ~] = fileparts(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}.fmrspace);
        [~, subrootb, ~] = fileparts(subroota);
        mat_file = sprintf('%s/%s.mat', foldera, subrootb);
        
        %If .mat file does not exist, create it
        if exist(mat_file, 'file') ~= 2
            mri = MRIread(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}.fmrspace);
            save(mat_file, 'mri');
            clear mri;
        end
        
        region = load(mat_file);
        fprintf('%s  Running ROI Analysis for ROI: %s\n', datestr(now), mat_file);diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
        %Create a 1D vector of voxel indices
        voxels = find(region.mri.vol >= 0.5); % Grab voxels with value > 0.5 (values below 1 represent blur from the image transformation)
        clear region;
        %Limit the number of voxels used to match the longitudinal group minimum for that ROI. Randomly select voxels to be thrown out.
        if isfield(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}, 'group_num_voxels')
            temp = randperm(length(voxels));
            voxels = sort(voxels(temp(1:data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}.group_num_voxels)));
            clear temp;
        end
        
        % Extract the ROI activation matrices.
        roi_matrix1 = matrix1(:,voxels);
        roi_matrix2 = matrix2(:,voxels);
        clear voxels;
        
        %Remove any voxels who's columns contain zeros
        inds = all(roi_matrix1)&all(roi_matrix2);
        roi_matrix1 = (roi_matrix1(:,inds));
        roi_matrix2 = (roi_matrix2(:,inds));
        
        %Store additional data for later use
        if ~isempty(datavar2)
            datavar2.apAnalyses(isess).apROI{kk}.activation.ave1 = ...
                mean(roi_matrix1, 2);
            datavar2.apAnalyses(isess).apROI{kk}.activation.ave2 = ...
                mean(roi_matrix2, 2);
            datavar2.apAnalyses(isess).apROI{kk}.activation.std1 = ...
                std(roi_matrix1, 0, 2);
            datavar2.apAnalyses(isess).apROI{kk}.activation.std2 = ...
                std(roi_matrix2, 0, 2);
            datavar2.apAnalyses(isess).apROI{kk}.activation.std_err1 = ...
                datavar2.apAnalyses(isess).apROI{kk}.activation.std1 / size(roi_matrix1, 2);
            datavar2.apAnalyses(isess).apROI{kk}.activation.std_err2 = ...
                datavar2.apAnalyses(isess).apROI{kk}.activation.std2 / size(roi_matrix2, 2);
        end
        
        % Perform the machine learning analysis
        results = bb_MVPA_Classifier_Test(roi_matrix1, roi_matrix2, 'classify');
        
        % Store results
        if ~isempty(datavar2)
            datavar2.apAnalyses(isess).apROI{kk}.results = results;
        end
        
        % Fill in the confusion matrix
        if ~isempty(conf_mat_loc)
            datavar(ii).mvpa.apAnalyses(isess).apROI{kk}.confusion_matrix.acc(conf_mat_loc(1), conf_mat_loc(2)) = results.cond1_mean_acc;
            datavar(ii).mvpa.apAnalyses(isess).apROI{kk}.confusion_matrix.acc(conf_mat_loc(2), conf_mat_loc(1)) = results.cond2_mean_acc;
            datavar(ii).mvpa.apAnalyses(isess).apROI{kk}.confusion_matrix.pval(conf_mat_loc(1), conf_mat_loc(2)) = results.pval_cond1;
            datavar(ii).mvpa.apAnalyses(isess).apROI{kk}.confusion_matrix.pval(conf_mat_loc(2), conf_mat_loc(1)) = results.pval_cond2;
        end
        clear results; clear roi_matrix1; clear roi_matrix2;
    end
end
