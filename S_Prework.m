function data = S_Prework(data, subNum, debug)

% Collect run information that will be used to analyze the fMRI data

if ~debug
    data.sub(subNum).study.stats.iv(2).task_conds = [12,50];
end
for rr=1:length(data.sub(subNum).study.runs)
    rrr = data.sub(subNum).study.runs(rr);
    tp_dur = data.sub(subNum).study.mristudy.fmri.fmr(rrr).TR/1000;
    iblock = 0;
    for kk=1:length(data.sub(subNum).run(rrr).comp.ivone)
        if ismember(data.sub(subNum).run(rrr).comp.ivone(kk), data.sub(subNum).study.stats.iv(2).task_conds)
            iblock = iblock+1;
            % Find beggining and end of the timecourse segment that represents each block.
            data.sub(subNum).run(rrr).stats.seg(iblock).start = data.sub(subNum).run(rrr).comp.start(kk-1)/tp_dur+1;
            data.sub(subNum).run(rrr).stats.seg(iblock).end = data.sub(subNum).run(rrr).comp.start(kk+1)/tp_dur+...
                data.sub(subNum).run(rrr).comp.durs(kk+1)/tp_dur;
        end
    end
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
    if isfield(data.sub(subNum).study.mristudy.fmri, 'fmr_fmrspaceMASK')
        input_file = (data.sub(subNum).study.mristudy.fmri.fmr_fmrspaceMASK);
        [subdirect, subroota, ~] = fileparts(data.sub(subNum).study.mristudy.fmri.fmr_fmrspaceMASK);
        data.sub(subNum).study.mristudy.fmri.mvpa.mask = sprintf('%s/%s.mat',subdirect, subroota);
    else
        input_file = sprintf('%s/%s/atlas/%s_fmr_fmrspace_ave_mask.nii', ...
            data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
        data.sub(subNum).study.mristudy.fmri.mvpa.mask = sprintf('%s/%s/atlas/%s_fmr_fmrspace_ave_mask.mat', ...
            data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
    end
    if exist(data.sub(subNum).study.mristudy.fmri.mvpa.mask, 'file') ~= 2
        mri = MRIread(input_file);
%         len = numel(mri.vol);
%         temp = sort(reshape(mri.vol, len,1));
%         cutoff = temp(round(len*.1))+(temp(round(len*.98))-temp(round(len*.1)))*1/4;
%         mri.vol(mri.vol<cutoff)=0;
%         %mri.vol(mri.vol<80)=0;
%         mri.vol(mri.vol~=0)=1;
        save(data.sub(subNum).study.mristudy.fmri.mvpa.mask, 'mri');
    end
end