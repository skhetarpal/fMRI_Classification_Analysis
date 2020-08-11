function [data, activation] = S_Data_Wrangling(study, data, subNum, debug)

%% Load the brain mask
mask = load(data.sub(subNum).study.mristudy.fmri.mvpa.mask);
data.sub(subNum).study.mristudy.fmri.mvpa.image_dims = size(mask.mri.vol);
image_size = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims;

%% Load fMRI data and calculate activation by run

for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    len = length(data.sub(subNum).run(rrr).stats.seg);

    %Load the run's fMRI data.
    clear mri; 
    mri = load(data.sub(subNum).study.mristudy.fmri.fmr(rrr).mcmat);
    fprintf('%s    %s\n', datestr(now), data.sub(subNum).study.mristudy.fmri.fmr(rrr).mcmat); if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
    
    %Initialize the activation matrices
    act.run(rrr).stim = zeros(image_size(1),image_size(2),image_size(3), len);
    act.run(rrr).null = zeros(image_size(1),image_size(2),image_size(3), len);
    
    %Load the timecourse for each voxel and calculate the activation
    for k=1:image_size(3)
        for j=1:image_size(2)
            for i=1:image_size(1)
                if mask.mri.vol(i,j,k) ~= 0
                    timecourse(:,1) = mri.mri.vol(i,j,k,:);
                    for seg = 1:len % Break timecourses into blocks.
                        segment = timecourse(data.sub(subNum).run(rrr).stats.seg(seg).start:data.sub(subNum).run(rrr).stats.seg(seg).end);
                        segment = zscore(segment);
                        if length(segment) == 25
                            tc_min = [4,5,6,24,25];
                            tc_max = [12,13,14,15,16];
                        elseif length(segment) == 30
                            tc_min = [9,10,11,29,30];
                            tc_max = [17,18,19,20,21];
end
                            act.run(rrr).stim(i,j,k,seg) = mean(segment(tc_max));
                            act.run(rrr).null(i,j,k,seg) = mean(segment(tc_min));
                    clear segment;
                    end
                    clear timecourse;
                end
            end
        end
    end
end

%Calculate the total number of blocks used in the analysis
total_len = 0;
for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    total_len = total_len + length(data.sub(subNum).run(rrr).stats.seg);
end

%Combine run data matrices into session data matrix
activation.stim = zeros(image_size(1),image_size(2),image_size(3),total_len);
activation.null = zeros(image_size(1),image_size(2),image_size(3),total_len);
counter = 1;
for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    run_length = size(act.run(rrr).stim,4);
    for k=1:image_size(3)
        for j=1:image_size(2)
            for i=1:image_size(1)
                activation.stim(i,j,k,counter:(counter+run_length-1)) = act.run(rrr).stim(i,j,k,:);
                activation.null(i,j,k,counter:(counter+run_length-1)) = act.run(rrr).null(i,j,k,:);
            end
        end
    end
    counter = counter + run_length;
end
clear act;

%% Reduce the size of the data matrices
activation.stim = single(activation.stim);
activation.null = single(activation.null);
    
