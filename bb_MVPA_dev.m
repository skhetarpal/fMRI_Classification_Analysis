function data = bb_MVPA_dev(study, data, subNum, html, upload, forprint)
 
%% Usage and Setup
fxname = mfilename;
if (nargin < 6)
    fprintf('\n');
    fprintf('Usage: data = %s(study, data, subNum, html, upload, forprint)\n', fxname);
    fprintf('       data = %s(study, data, subNum, 1,    1,      1)\n', fxname);
    fprintf('study contains important variables for the study.\n');
    fprintf('study should be contained within a study-level parameter file.\n');
    fprintf('data contains important variables for the subject.\n');
    fprintf('data should be contained within a study-level parameter file.\n');
    fprintf('subNum is the subject number.\n');
    fprintf('html is whether (1) or not (0) to create the html.\n');
    fprintf('upload is either 1 indicating upload all data to server or 0 indicating do not upload\n');
    fprintf('forprint is whether (1) or not (0) to format for publication.\n');
    return;
end

msg = nargoutchk(1, 1, nargout);
if ~isempty(msg)
    fprintf('%s: %s', fxname,msg);
    return;
end

addpath /autofs/cluster/vincent/programs/matlab/bb_toolbox_130701;
addpath /autofs/cluster/vincent/programs/matlab/suraj;
close all;
tic; % start time
chdir(fullfile(study.info.studydir, data.sub(subNum).info.name)); % Enter Subject Directory
% Create MVPA Directories
data.sub(subNum).qc.MVPAaaqcdir = sprintf('%s/%s/%s_MVPAaa_qc', data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
if (exist(data.sub(subNum).qc.MVPAaaqcdir, 'dir') ~= 7)
    mkdir(data.sub(subNum).qc.MVPAaaqcdir);
else
    cmd = sprintf('!rm %s/*',data.sub(subNum).qc.MVPAaaqcdir);
    eval(cmd);
end


% %Load ROIs from bb_MVPA_Group and transform them to native space.
% fprintf('Running bb_MVPA_ROI_Prep on %s\n', data.sub(subNum).info.name);
% data = bb_MVPA_ROI_Prep(study, data, subNum, 0);

% %Create ROI overlays.
% fprintf('Running bb_MVPA_Create_ROI_Overlays on %s\n', data.sub(subNum).info.name);
% data = bb_MVPA_Create_ROI_Overlays(study, data, subNum);

% %Load the brain mask, transform it to native space, and create the mask overlays.
% fprintf('Running bb_MVPA_Mask_Prep on %s\n', data.sub(subNum).info.name);
% data = bb_MVPA_Mask_Prep(data, subNum);

%Create the regression matrix, and create mat files from each run's fMRI data.
fprintf('Running bb_MVPA_Prework on %s\n', data.sub(subNum).info.name);
data = bb_MVPA_Prework(data, subNum);

%Load each run's fMRI data and calculate beta activations as well as task effect activations (stimulus vs. null).
fprintf('Running bb_MVPA_Data_Wrangling on %s\n', data.sub(subNum).info.name);
[data, activation] = bb_MVPA_Data_Wrangling(study, data, subNum);

%Check that consecutive blocls are independent from each other.
fprintf('Running bb_MVPA_Independence_Test on %s\n', data.sub(subNum).info.name);
data = bb_MVPA_Independence_Test(data, subNum, activation.session.betas);

%Loop through each level of bv, iv, dv, and mv.  Perform searchlight analyses as well as ROI analyses.
fprintf('Running bb_MVPA_Analysis on %s\n', data.sub(subNum).info.name);
data = bb_MVPA_Analysis(study, data, subNum, activation);

datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
save(datafile, 'data');

%%%%%%%%%%%%%%%
% Generate HTML
%%%%%%%%%%%%%%%
if (html)
    fprintf('\nGenerating HTML Markup.\n');
    if (upload)
        if ~exist(study.export.serverdir, 'dir')
            mkdir(study.export.serverdir)
        end
        chdir(study.export.serverdir);
    else
        chdir(data.sub(subNum).qc.MVPAaaqcdir);
    end
    
    mainhtmlfilename    = sprintf('%s_%s_mvpa_report.html', study.info.id, data.sub(subNum).info.name);
    htmltitle           = sprintf('%s\n%s MVPA Report', study.info.id, data.sub(subNum).info.name);
    fidhtml             = bb_open_html(study, mainhtmlfilename, htmltitle);

    % HTML: Generate Table of Contents
    fidhtml = bb_HTMLSPMFunctionalAnalysisToC(study, data, subNum, upload, fidhtml, forprint);

    % HTML: Study Information
    fidhtml = bb_HTMLSPMFunctionalAnalysisPrintStudyInformation(study, data, subNum, fidhtml);
    fclose(fidhtml);
    
    % HTML: Quality Control Overview
    data = bb_HTMLSPMFunctionalAnalysisQCOverview(study, data, subNum, mainhtmlfilename, 'MVPA', upload, forprint);

    %%%
    % HTML: MVPA Contrasts
    %%%
    data = bb_HTMLMVPAContrasts(study, data, subNum, mainhtmlfilename, forprint);
    
    % Close HTML
    fidhtml = fopen(mainhtmlfilename, 'a');
    x = toc; % elapsed program time
    bb_close_html(mfilename, fidhtml, forprint, x);
end

datafile = sprintf('%s/%s/data_%s_%s_bb_MVPA_fmri.mat', study.info.studydir, data.sub(subNum).info.name, study.info.id, data.sub(subNum).info.name);
save(datafile, 'data');

chdir(study.info.studydir);
fprintf('\n%s Complete for %s.\n', fxname, data.sub(subNum).info.name);