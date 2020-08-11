function datavar = bb_MVPA_Searchlight_Atlas_Transformation(data, subNum, datavar)

for aa=1:length(data.sub(subNum).study.mristudy.fmratltargs)
    aaa = data.sub(subNum).study.mristudy.fmratltargs(aa);
    target = data.sub(subNum).study.mristudy.atl(aaa).target;
    
    if (data.sub(subNum).study.mristudy.atl(aaa).method == 1)
        % FSL Registration
        datavar.searchlight.atl(aaa).fname = ...
            sprintf('%s_atl%d', datavar.searchlight_file, aaa);
        cmd = sprintf('!flirt -ref %s -in %s.nii -out %s.nii -nosearch -applyxfm -init %s -interp trilinear', ...
            target, ...
            datavar.searchlight_file, ...
            datavar.searchlight.atl(aaa).fname, ...
            data.sub(subNum).study.mristudy.mat(aaa).fmr2target);
        eval(cmd);
        
    elseif (data.sub(subNum).study.mristudy.atl(aaa).method == 2)
        if isfield(data.sub(subNum).study.mristudy, 'com')
            if exist(data.sub(subNum).study.mristudy.com(aaa).fmr2target, 'file')
                % JIP Registration
                datavar.searchlight.atl(aaa).fname = ...
                    sprintf('%s_atl%d', datavar.searchlight_file, aaa);
                cmd = sprintf('!%s/jip %s/AlignNii.com %s.nii %s %s.nii', ...
                    data.sub(subNum).study.mristudy.info.jippath, ...
                    data.sub(subNum).study.mristudy.info.jippath, ...
                    datavar.searchlight_file, ...
                    data.sub(subNum).study.mristudy.com(aaa).fmr2target, ...
                    datavar.searchlight.atl(aaa).fname);
                eval(cmd);
            else
                fprintf('\n %s: no JIP\n', data.sub(subNum).info.name);
            end
        else
            fprintf('\n %s: no JIP\n', data.sub(subNum).info.name);
        end
    end
end
end