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
    session(1).name = 'fMRI_baby4_140713';
    session(2).name = 'fMRI_baby4_140727';
    session(3).name = 'fMRI_baby4_141011';
    session(4).name = 'fMRI_baby4_141109';
    session(5).name = 'fMRI_baby4_141123';
    session(6).name = 'fMRI_baby4_141224';
    session(7).name = 'fMRI_baby4_150222';
    session(8).name = 'fMRI_baby4_150524';
    session(9).name = 'fMRI_baby4_150705';
    study_(1).name = 'Study000012';
    study_(2).name = 'Study000012';
    study_(3).name = 'Study000012';
    study_(4).name = 'Study000012';
    study_(5).name = 'Study000012';
    study_(6).name = 'Study000012';
    study_(7).name = 'Study000012';
    study_(8).name = 'Study000012';
    study_(9).name = 'Study000012';
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
        else
            cmd = sprintf('!rm %s/*',data.sub(subNum).qc.MVPAaaqcdir);
            eval(cmd);
        end
    end
    
    x(ss).age = data.sub(subNum).info.age;
    for ap=1:3%length(study.part(data.sub(subNum).info.part).newapAnalysis)
        isess = study.part(data.sub(subNum).info.part).newapAnalysis(ap);
        for kk=1:length(data.sub(subNum).study.part(data.sub(subNum).info.part).apAnalyses(isess).apROI)
            x(ss).apAnalyses(isess).apROI{kk}.acc = data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.results.acc;
            x(ss).apAnalyses(isess).apROI{kk}.pval = data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.results.pval;
x(ss).apAnalyses(isess).apROI{kk}.perm_upper = data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.results.perm_upper;
x(ss).apAnalyses(isess).apROI{kk}.num_obs = data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.num_obs;
x(ss).apAnalyses(isess).apROI{kk}.num_vox = data.sub(subNum).study.stats.iv(2).mvpa.apAnalyses(isess).apROI{kk}.num_vox;
        end
    end
end

for ss=1:length(session)
    age(ss) = x(ss).age;
    acc(ss) = x(ss).apAnalyses(3).apROI{1}.acc;
    pval(ss) = x(ss).apAnalyses(3).apROI{1}.pval;
perm_upper(ss) = x(ss).apAnalyses(3).apROI{1}.perm_upper;
end 
plot(age,acc,'rx',age,perm_upper,'gx'); ylim([.5,1.1]);
