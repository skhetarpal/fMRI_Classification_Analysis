function data = bb_MVPA_Prework(data, subNum)

%% Make MION HRF
%This produces a vector that represents a typical Hemodynamic Response
if ~debug
    hrf = bb_MakeHRF(1/16, 513/16,  1.5,  4.5, 13.5,   -0.18,    0.33,    0.67);
    spacing = 32;
    h.response = hrf.response(1:spacing:end); %Rescale it to 0.5 hertz
end

%% Collect run information that will be used to analyze the fMRI data

for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    tp_dur = data.sub(subNum).study.mristudy.fmri.fmr(rrr).TR/1000;
    ntimepoints = sum(data.sub(subNum).run(rrr).comp.durs)/tp_dur;
    iblock = 0;
    for kk=1:length(data.sub(subNum).run(rrr).comp.ivone)
        if ismember(data.sub(subNum).run(rrr).comp.ivone(kk), data.sub(subNum).study.stats.iv(1).plotord)
            iblock = iblock+1;
            data.sub(subNum).run(rrr).comp.cond_sequence(iblock) = data.sub(subNum).run(rrr).comp.ivone(kk);
            
            % Find beggining and end of the timecourse segment that represents each block.
            data.sub(subNum).run(rrr).stats.seg(iblock).start = data.sub(subNum).run(rrr).comp.start(kk-1)/tp_dur+1;
            data.sub(subNum).run(rrr).stats.seg(iblock).end = data.sub(subNum).run(rrr).comp.start(kk+1)/tp_dur+...
                data.sub(subNum).run(rrr).comp.durs(kk+1)/tp_dur;
            
            % Setup the hrf convolved matrix for the regression analysis
            on = data.sub(subNum).run(rrr).comp.start(kk)/tp_dur+1;
            off = on + data.sub(subNum).run(rrr).comp.durs(kk)/tp_dur-1;
            data.sub(subNum).run(rrr).stats.cm(iblock).vector = zeros(ntimepoints, 1);
            data.sub(subNum).run(rrr).stats.cm(iblock).vector(on:off) = 1; 
            data.sub(subNum).run(rrr).stats.cm(iblock).vector = conv(data.sub(subNum).run(rrr).stats.cm(iblock).vector, h.response);
            data.sub(subNum).run(rrr).stats.cm(iblock).vector = data.sub(subNum).run(rrr).stats.cm(iblock).vector(1:ntimepoints);
            % Remove the variance. (Removing the mean would change the x-intercepts, which is undesirous)
            data.sub(subNum).run(rrr).stats.cm(iblock).vector = data.sub(subNum).run(rrr).stats.cm(iblock).vector/std(data.sub(subNum).run(rrr).stats.cm(iblock).vector);
        end
    end
    % Create the run's regression matrix "XX"
    XX = [];
    %Add expected response columns for each condition present during the run
    for iblock=1:length(data.sub(subNum).run(rrr).stats.cm)
        XX = [XX data.sub(subNum).run(rrr).stats.cm(iblock).vector];
    end
    XX = [XX ones(ntimepoints,1)]; %Add in the Baseline Nuissance Regressor
    XX = [XX ((1:ntimepoints)' - (ntimepoints+1)/2)]; %Add in the Trend Nuissance Regressor
    %Add in a Motion Nuissance Regressor
    t0 = data.sub(subNum).study.mristudy.fmri.fmr(rr).nvol*(rr-1)+1; %Start of run
    t1 = t0+data.sub(subNum).study.mristudy.fmri.fmr(rr).nvol-1; %End of the run
    for zz=1:size(data.sub(subNum).qc.motionRegressors,2)
        if sum(abs(data.sub(subNum).qc.motionRegressors(t0:t1,zz))) > 0
            XX = [XX data.sub(subNum).qc.motionRegressors(t0:t1, zz)];
        end
    end
    data.sub(subNum).run(rrr).stats.XX_inv = pinv(XX);
%     data.sub(subNum).run(rrr).stats.pinv_XTX_XT = pinv(transpose(XX)*XX)*transpose(XX);
end

%Create a master vector containing the combined ivone sequence from all runs included in the analysis.
data.sub(subNum).cond_sequence=[];
for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    data.sub(subNum).cond_sequence = [data.sub(subNum).cond_sequence, data.sub(subNum).run(rrr).comp.cond_sequence];
end

%% Check if each run's fMRI .mat file exists. If not, create it.
if ~debug
    for rr=1:length(data.sub(subNum).study.runs)
        rrr = data.sub(subNum).study.runs(rr);
        [subdirect, subroota, ~] = fileparts(data.sub(subNum).study.mristudy.fmri.fmr(rrr).mc);
        [~, subrootb, ~] = fileparts(subroota);
        data.sub(subNum).study.mristudy.fmri.fmr(rrr).mcmat = sprintf('%s/%s.mat', subdirect, subrootb);
        if exist(data.sub(subNum).study.mristudy.fmri.fmr(rrr).mcmat, 'file') ~= 2
            mri = MRIread(data.sub(subNum).study.mristudy.fmri.fmr(rrr).mc);
            save([subdirect '/' subrootb '.mat'], 'mri')
        end
    end
end

%% If the mask does not exist, create it
if ~debug
    data.sub(subNum).study.mristudy.fmri.mvpa.mask = sprintf('%s/%s/atlas/%s_mvpa_mask_iv2task.nii', ...
        data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
    if exist(data.sub(subNum).study.mristudy.fmri.mvpa.mask, 'file') ~= 2
        input_file = sprintf('%s/%s/spm_results/IV2_Task/%s_spmT_Task_UNC0p001_posneg.nii', ...
            data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
        mri = MRIread(input_file);
        mri.vol(mri.vol~=0)=1;
        save(data.sub(subNum).study.mristudy.fmri.mvpa.mask, 'mri');
    end
end