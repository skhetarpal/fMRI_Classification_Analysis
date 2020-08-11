function [data, activation] = bb_MVPA_Data_Wrangling(study, data, subNum)

%% Load the brain mask
mask = load(data.sub(subNum).study.mristudy.fmri.mvpa.mask);
data.sub(subNum).study.mristudy.fmri.mvpa.image_dims = size(mask.mri.vol);
image_size = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims;
num_voxels = numel(mask.mri.vol);

%% Load fMRI data and calculate activation by run

for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    
    %Load the run's fMRI data.
    clear mri; 
    mri = load(data.sub(subNum).study.mristudy.fmri.fmr(rrr).mcmat);
    fprintf('%s    %s\n', datestr(now), data.sub(subNum).study.mristudy.fmri.fmr(rrr).mcmat); diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
    
%    if rr==1
%        data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(1) = size(mri.mri.vol,1); 
%        data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(2) = size(mri.mri.vol,2); 
%        data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(3) = size(mri.mri.vol,3);
%        num_voxels = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(1) * data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(2) * data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(3);
%    end
    
    %Initialize the activation matrices
    act.run(rrr).stim_act = zeros(length(data.sub(subNum).run(rrr).stats.seg), num_voxels);
    act.run(rrr).null_act = zeros(length(data.sub(subNum).run(rrr).stats.seg), num_voxels);
    act.run(rrr).betas = zeros(size(data.sub(subNum).run(rrr).stats.XX_inv,1), num_voxels);
    
    %Load the timecourse for each voxel and calculate the activation
    counter = 1;
    for k=1:data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(3)
        for j=1:data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(2)
            for i=1:data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(1)
                if mask.mri.vol(i,j,k) ~= 0
                    timecourse(:,1) = mri.mri.vol(i,j,k,:);
                    voxel_betas = data.sub(subNum).run(rrr).stats.XX_inv * timecourse;
%                     voxel_betas = data.sub(subNum).run(rrr).stats.pinv_XTX_XT * timecourse;
                    act.run(rrr).betas(:,counter) = voxel_betas;
                    for seg = 1:length(data.sub(subNum).run(rrr).stats.seg) % Break timecourses into blocks.
                        segment = timecourse(data.sub(subNum).run(rrr).stats.seg(seg).start:data.sub(subNum).run(rrr).stats.seg(seg).end);
                        segment = zscore(segment);
                        act.run(rrr).stim_act(seg,counter) = mean(segment(study.mristudy.avgtc.tc_max));
                        act.run(rrr).null_act(seg,counter) = mean(segment(study.mristudy.avgtc.tc_min));
                    end
                    clear segment;
                    clear timecourse;
                end
                counter = counter + 1;
            end
        end
    end
end

%% Combine run data matrices into session data matrix

activation.session.betas = zeros(length(data.sub(subNum).cond_sequence), num_voxels);
activation.session.stim_act = zeros(length(data.sub(subNum).cond_sequence), num_voxels);
activation.session.null_act = zeros(length(data.sub(subNum).cond_sequence), num_voxels);
counter = 1;
for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    run_length = length(data.sub(subNum).run(rrr).comp.cond_sequence);
    activation.session.betas(counter:(counter+run_length-1),:) = act.run(rrr).betas(1:run_length, :);
    activation.session.stim_act(counter:(counter+run_length-1),:) = act.run(rrr).stim_act;
    activation.session.null_act(counter:(counter+run_length-1),:) = act.run(rrr).null_act;
    counter = counter + run_length;
end
clear act;

%% Reduce the size of the data matrices
activation.session.betas = single(activation.session.betas);
activation.session.stim_act = single(activation.session.stim_act);
activation.session.null_act = single(activation.session.null_act);
    
