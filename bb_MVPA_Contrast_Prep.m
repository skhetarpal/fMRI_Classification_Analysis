function datavar = bb_MVPA_Contrast_Prep(data, subNum, datavar, ii)

for icontrast=1:length(datavar(ii).contrast)
    if isfield(datavar(ii).contrast(icontrast).tcon, 'suffTrials') && datavar(ii).contrast(icontrast).tcon.suffTrials
        datavar(ii).contrast(icontrast).mvpa.suffTrials = 1;
        %Confirm that there are enough trails to do an analysis
        for cc=1:length(datavar(ii).contrast(icontrast).tcon.convec)
            if datavar(ii).contrast(icontrast).tcon.convec(cc) ~= 0
                %Confirm that each level in the contrast has a sufficient number of blocks.
                for tt = 1:length(datavar(ii).level(cc).types)
                    temp(tt,:) = datavar(ii).level(cc).types(tt) == data.sub(subNum).cond_sequence;
                end
                if exist('temp') && sum(sum(temp)) < 4
                    datavar(ii).contrast(icontrast).mvpa.suffTrials = 0;
                end
            end
        end
    else
        datavar(ii).contrast(icontrast).mvpa.suffTrials = 0;
    end
    
    if datavar(ii).contrast(icontrast).mvpa.suffTrials
        datavar(ii).contrast(icontrast).mvpa.name = datavar(ii).contrast(icontrast).tcon.name;
        datavar(ii).contrast(icontrast).mvpa.convec = datavar(ii).contrast(icontrast).tcon.convec;
        
        
        datavar(ii).contrast(icontrast).mvpa.cond1_levels = [];
        datavar(ii).contrast(icontrast).mvpa.cond2_levels = [];
        datavar(ii).contrast(icontrast).mvpa.cond1_name = [];
        datavar(ii).contrast(icontrast).mvpa.cond2_name = [];
        for cc=1:length(datavar(ii).contrast(icontrast).mvpa.convec)
            if datavar(ii).contrast(icontrast).mvpa.convec(cc) > 0
                datavar(ii).contrast(icontrast).mvpa.cond1_levels = [datavar(ii).contrast(icontrast).mvpa.cond1_levels cc];
                if isempty(datavar(ii).contrast(icontrast).mvpa.cond1_name)
                    datavar(ii).contrast(icontrast).mvpa.cond1_name = datavar(ii).level(cc).name;
                else
                    datavar(ii).contrast(icontrast).mvpa.cond1_name = strcat(datavar(ii).contrast(icontrast).mvpa.cond1_name, '_and_', datavar(ii).level(cc).name);
                end
            elseif datavar(ii).contrast(icontrast).mvpa.convec(cc) < 0
                datavar(ii).contrast(icontrast).mvpa.cond2_levels = [datavar(ii).contrast(icontrast).mvpa.cond2_levels cc];
                if isempty(datavar(ii).contrast(icontrast).mvpa.cond2_name)
                    datavar(ii).contrast(icontrast).mvpa.cond2_name = datavar(ii).level(cc).name;
                else
                    datavar(ii).contrast(icontrast).mvpa.cond2_name = strcat(datavar(ii).contrast(icontrast).mvpa.cond2_name, '_and_', datavar(ii).level(cc).name);
                end
           end
        end
        datavar(ii).contrast(icontrast).mvpa.cond1_blocks = ismember(data.sub(subNum).cond_sequence, datavar(ii).contrast(icontrast).mvpa.cond1_levels);
        datavar(ii).contrast(icontrast).mvpa.cond2_blocks = ismember(data.sub(subNum).cond_sequence, datavar(ii).contrast(icontrast).mvpa.cond2_levels);
        datavar(ii).contrast(icontrast).mvpa.long_name = strcat(datavar(ii).contrast(icontrast).mvpa.cond1_name, '_vs_', datavar(ii).contrast(icontrast).mvpa.cond2_name);
    end
end
end