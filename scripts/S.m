dbstop if error;
debug = 0;
if ~debug
    diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
    addpath /autofs/cluster/vincent/programs/matlab/bb_toolbox_130701;
    addpath /autofs/cluster/vincent/programs/matlab/suraj/7_10_16;
end
html=1;
upload = 1;
forprint = 1;

if debug
    session(1).name = 'test_name';
    study_(1).name = 'test_study';
else
    %session(1).name = 'fMRI_baby4_140713';
    %session(2).name = 'fMRI_baby4_140727';
    %session(3).name = 'fMRI_baby4_141011';
    %session(4).name = 'fMRI_baby4_141109';
    %session(5).name = 'fMRI_baby4_141123';
    %session(6).name = 'fMRI_baby4_141224';
    %session(7).name = 'fMRI_baby4_150222';
    %session(8).name = 'fMRI_baby4_150524';
    session(1).name = 'fMRI_baby4_150705';
session(2).name = 'fMRI_baby4_140713';
    %session(10).name = 'fMRI_baby5_130706';
    %session(1).name = 'fMRI_baby5_141011';
    %session(2).name = 'fMRI_baby5_141109';
    %session(3).name = 'fMRI_baby5_141123';
    %session(4).name = 'fMRI_baby5_141224';
    %session(5).name = 'fMRI_baby5_150125';
    %session(6).name = 'fMRI_baby5_150222';
    %session(7).name = 'fMRI_baby5_150524';
    %session(8).name = 'fMRI_baby5_150628';
    study_(1).name = 'Study000012';
    study_(2).name = 'Study000012';
    study_(3).name = 'Study000012';
    study_(4).name = 'Study000012';
    study_(5).name = 'Study000012';
    study_(6).name = 'Study000012';
    study_(7).name = 'Study000012';
    study_(8).name = 'Study000012';
    study_(9).name = 'Study000012';
%     study_(10).name = 'Study000012';
%     study_(11).name = 'Study000012';
%     study_(12).name = 'Study000012';
%     study_(13).name = 'Study000012';
%     study_(14).name = 'Study000012';
%     study_(15).name = 'Study000012';
%     study_(16).name = 'Study000012';
%     study_(17).name = 'Study000012';
%     study_(18).name = 'Study000012';
end

%Run Analysis of Sessions
for ss=1:length(session)
    if debug
        [study, data, subNum] = Test_Environment();
    else
        clear study; clear data;
        study.info.id = study_(ss).name;
        run(sprintf('/autofs/cluster/vincent/livingstone/%s/%s_pp_params.m', study.info.id, study.info.id));
        run(sprintf('/autofs/cluster/vincent/livingstone/%s/%s_aa_params.m', study.info.id, study.info.id));
        run(sprintf('/autofs/cluster/vincent/livingstone/%s/%s_Data/%s_%s_pp_params.m', ...
            study.info.id, session(ss).name, study.info.id, session(ss).name));
        run(sprintf('/autofs/cluster/vincent/livingstone/%s/%s_Data/%s_%s_aa_params.m', ...
            study.info.id, session(ss).name, study.info.id, session(ss).name));
        datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
        if exist(datafile) == 2
            fprintf('Loading MVPA file for session %s\n', session(ss).name); diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
            load(datafile);
        else
            input_file = sprintf('/autofs/cluster/vincent/livingstone/%s/%s/data_%s_%s_bb_FunctionalAnalysis_fmri.mat', ...
                study.info.id, session(ss).name, study.info.id, session(ss).name);
            if exist(input_file) == 2
                fprintf('Loading Functional Analysis file for session %s\n', session(ss).name); diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
                load(input_file);
            end
        end
        data.sub(subNum).qc.MVPAaaqcdir = sprintf('%s/%s/%s_MVPAaa_qc', data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
        if (exist(data.sub(subNum).qc.MVPAaaqcdir, 'dir') ~= 7)
            mkdir(data.sub(subNum).qc.MVPAaaqcdir);
        end
    end
    
    fprintf('%s  Running bb_MVPA_Prework on %s\n', datestr(now), data.sub(subNum).info.name); if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
    data = S_Prework(data, subNum, debug);
    fprintf('%s  Running bb_MVPA_Data_Wrangling on %s\n', datestr(now), data.sub(subNum).info.name); if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
    [data, activation] = S_Data_Wrangling(study, data, subNum, debug);
    fprintf('%s  Running S_Independence_Test on %s\n', datestr(now), data.sub(subNum).info.name); if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
    data = S_Independence_Test(data, subNum, activation);
    fprintf('%s  Running bb_MVPA_ROI_Analysis on %s\n', datestr(now), data.sub(subNum).info.name);
    data = S_ROI_Analysis(study, data, subNum, activation, debug);
if ss == 1
    fprintf('%s  Running bb_MVPA_Searchlight_Analysis on %s\n', datestr(now), data.sub(subNum).info.name);
    data = S_Searchlight_Analysis(data, subNum, activation, debug);
    if debug
        datafile = sprintf('%s/data_%s_%s_bb_MVPA_fmri.mat', data.sub(subNum).qc.MVPAaaqcdir, study.info.id, data.sub(subNum).info.name);
    else
        datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
    end
    save(datafile, 'data');
    clear data;
end
end
