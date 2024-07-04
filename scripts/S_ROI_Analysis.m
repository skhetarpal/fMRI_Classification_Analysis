function data = S_ROI_Analysis(study, data, subNum, activation, debug)

for ap=1:3%length(study.part(data.sub(subNum).info.part).newapAnalysis)
    isess = study.part(data.sub(subNum).info.part).newapAnalysis(ap);
    for kk=1:length(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI)
        if debug
            mat_file = data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}.fmrspace;
        else
            %Determine the name of the .mat file
            [foldera, subroota, ~] = fileparts(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}.fmrspace);
            [~, subrootb, ~] = fileparts(subroota);
            mat_file = sprintf('%s/%s.mat', foldera, subrootb);
            
            %If .mat file does not exist, create it
            if exist(mat_file, 'file') ~= 2
                mri = MRIread(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}.fmrspace);
                save(mat_file, 'mri');
                clear mri;
            end
        end
        
        region = load(mat_file);
        fprintf('%s  Running ROI Analysis for ROI: %s\n', datestr(now), mat_file);if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
        %Create a matrix of logicals for voxels in the ROI
        voxels = region.mri.vol >= 0.5; % Grab voxels with value > 0.5 (values below 1 represent blur from the image transformation)
        clear region;
        num_voxels = nnz(voxels);
        len = size(activation.stim,4);
        roi_stim = zeros(len, num_voxels);
        roi_null = zeros(len, num_voxels);
        image_size = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims;
        counter = 0;
        for k=1:image_size(3)
            for j=1:image_size(2)
                for i=1:image_size(1)
                    if voxels(i,j,k)
                        if all(activation.stim(i,j,k,:)) && all(activation.null(i,j,k,:))
                            counter = counter + 1;
                            roi_stim(:,counter) = activation.stim(i,j,k,:);
                            roi_null(:,counter) = activation.null(i,j,k,:);
                        end
                    end
                end
            end
        end
        roi_stim = roi_stim(:,roi_stim(1,:)~=0);
        roi_null = roi_null(:,roi_stim(1,:)~=0);
        clear voxels;
        %Limit the number of voxels used to match the longitudinal group minimum for that ROI. Randomly select voxels to be thrown out.
        if isfield(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}, 'group_num_voxels')
            temp = randperm(counter);
            temp = sort(temp(1:data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI{kk}.group_num_voxels));
            roi_stim = roi_stim(:,temp);
            roi_stim = roi_stim(:,temp);
            clear temp;
        end
        
        %Store additional data for later use
	data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.num_obs = len * 2;
	data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.num_vox = size(roi_stim, 2);
        data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.ave_stim = mean(roi_stim, 2);
        data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.activation.ave_null = mean(roi_null, 2);
        data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.activation.stim_std = std(roi_stim, 0, 2);
        data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.activation.null_std = std(roi_null, 0, 2);
        data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.activation.stim_std_err = ...
            data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.activation.stim_std / size(roi_stim, 2);
        data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.activation.null_std_err = ...
            data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.activation.null_std / size(roi_null, 2);
        
        % Perform the machine learning analysis
        data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.results = ...
            bb_MVPA_Classifier_Test(roi_stim, roi_null, 'classify');
        clear roi_stim; clear roi_null;
        fprintf('%s  Accuracy for ROI %s is %f\n', datestr(now), mat_file, ...
            data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.results.acc);if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
    end
end
