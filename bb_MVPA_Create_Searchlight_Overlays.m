function datavar = bb_MVPA_Create_Searchlight_Overlays(data, subNum, datavar, threshp)

threshstat = 'anal';

output = sprintf('%s_overlay_%4.2f_%4.2f', datavar.searchlight.atl(data.sub(subNum).study.mristudy.fmratldisp).fname,...
    threshp, datavar.anal.threshmaxv);

cmd = sprintf('!overlay 0 1 %s %d %d %s %f %f %s;', ...
    data.sub(subNum).study.mristudy.atl(data.sub(subNum).study.mristudy.fmratldisp).fmr_on_target, ...
    data.sub(subNum).study.mristudy.fmri.thresholds(1), ...
    data.sub(subNum).study.mristudy.fmri.thresholds(2), ...
    datavar.searchlight.atl(data.sub(subNum).study.mristudy.fmratldisp).fname, ...
    threshp, ...
    datavar.anal.threshmaxv, ...
    output);
fprintf('%s\n', cmd);
eval(cmd);

% Saggital Left
cmd = sprintf('datavar.searchlight.fig.%s.lsag = ''%s_Lsag.gif'';', threshstat, output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!slicer %s.nii.gz -t -s 1 -x 0.10 sl10 -x 0.15 sl15 -x 0.20 sl20 -x 0.25 sl25 -x 0.30 sl30 -x 0.35 sl35 -x 0.40 sl40 -x 0.45 sl45 -x 0.50 sl50;', output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!convert -colors 100 +append sl10 sl15 sl20 sl25 sl30 sl35 sl40 sl45 sl50 %s_Lsag.gif;', output);
fprintf('%s\n', cmd);
eval(cmd);
!rm -f sl10 sl15 sl20 sl25 sl30 sl35 sl40 sl45 s150
% cmd = sprintf('movefile(datavar.searchlight.fig.%s.lsag, data.sub(subNum).qc.MVPAaaqcdir);', threshstat);
% fprintf('%s\n', cmd);
% eval(cmd);

% Saggital Right
cmd = sprintf('datavar.searchlight.fig.%s.rsag = ''%s_Rsag.gif'';', threshstat, output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!slicer %s.nii.gz -t -s 1 -x 0.55 sl55 -x 0.60 sl60 -x 0.65 sl65 -x 0.70 sl70 -x 0.75 sl75 -x 0.80 sl80 -x 0.85 sl85 -x 0.90 sl90;', output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!convert -colors 100 +append sl55 sl60 sl65 sl70 sl75 sl80 sl85 sl90 %s_Rsag.gif;', output);
fprintf('%s\n', cmd);
eval(cmd);

!rm -f sl55 sl60 sl65 sl70 sl75 sl80 sl85 sl90
% cmd = sprintf('movefile(datavar.searchlight.fig.%s.rsag, data.sub(subNum).qc.MVPAaaqcdir);', threshstat);
% fprintf('%s\n', cmd);
% eval(cmd);

% Transverse Ventral
cmd = sprintf('datavar.searchlight.fig.%s.vtra = ''%s_Vtra.gif'';', threshstat, output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!slicer %s.nii.gz -t -s 1 -z 0.15 sl15 -z 0.20 sl20 -z 0.25 sl25 -z 0.30 sl30 -z 0.35 sl35 -z 0.40 sl40 -z 0.45 sl45 -z 0.50 sl50;', output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!convert -colors 100 +append sl15 sl20 sl25 sl30 sl35 sl40 sl45 sl50 %s_Vtra.gif;', output);
fprintf('%s\n', cmd);
eval(cmd);

!rm -f sl15 sl20 sl25 sl30 sl35 sl40 sl45 sl50
% cmd = sprintf('movefile(datavar.searchlight.fig.%s.vtra, data.sub(subNum).qc.MVPAaaqcdir);', threshstat);
% fprintf('%s\n', cmd);
% eval(cmd);

% Transverse Dorsal
cmd = sprintf('datavar.searchlight.fig.%s.dtra = ''%s_Dtra.gif'';', threshstat, output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!slicer %s.nii.gz -t -s 1 -z 0.55 sl55 -z 0.60 sl60 -z 0.65 sl65 -z 0.70 sl70 -z 0.75 sl75 -z 0.80 sl80 -z 0.85 sl85;', output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!convert -colors 100 +append sl55 sl60 sl65 sl70 sl75 sl80 sl85 %s_Dtra.gif;', output);
fprintf('%s\n', cmd);
eval(cmd);

!rm -f sl55 sl60 sl65 sl70 sl75 sl80 sl85
% cmd = sprintf('movefile(datavar.searchlight.fig.%s.dtra, data.sub(subNum).qc.MVPAaaqcdir);', threshstat);
% fprintf('%s\n', cmd);
% eval(cmd);

% Coronal Posterior
cmd = sprintf('datavar.searchlight.fig.%s.pcor = ''%s_Pcor.gif'';', threshstat, output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!slicer %s.nii.gz -t -s 1 -y 0.10 sl10 -y 0.15 sl15 -y 0.20 sl20 -y 0.25 sl25 -y 0.30 sl30 -y 0.35 sl35 -y 0.40 sl40 -y 0.45 sl45 -y 0.50 sl50;', output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!convert -colors 100 +append sl10 sl15 sl20 sl25 sl30 sl35 sl40 sl45 sl50 %s_Pcor.gif;', output);
fprintf('%s\n', cmd);
eval(cmd);

!rm -f sl10 sl15 sl20 sl25 sl30 sl35 sl40 sl45 sl50
% cmd = sprintf('movefile(datavar.searchlight.fig.%s.pcor, data.sub(subNum).qc.MVPAaaqcdir);', threshstat);
% fprintf('%s\n', cmd);
% eval(cmd);

% Coronal Anterior
cmd = sprintf('datavar.searchlight.fig.%s.acor = ''%s_Acor.gif'';', threshstat, output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!slicer %s.nii.gz -t -s 1 -y 0.55 sl55 -y 0.60 sl60 -y 0.65 sl65 -y 0.70 sl70 -y 0.75 sl75 -y 0.80 sl80 -y 0.85 sl85;', output);
fprintf('%s\n', cmd);
eval(cmd);

cmd = sprintf('!convert -colors 100 +append sl55 sl60 sl65 sl70 sl75 sl80 sl85 %s_Acor.gif;', output);
fprintf('%s\n', cmd);
eval(cmd);

!rm -f sl55 sl60 sl65 sl70 sl75 sl80 sl85
% cmd = sprintf('movefile(datavar.searchlight.fig.%s.acor, data.sub(subNum).qc.MVPAaaqcdir);', threshstat);
% fprintf('%s\n', cmd);
% eval(cmd);
