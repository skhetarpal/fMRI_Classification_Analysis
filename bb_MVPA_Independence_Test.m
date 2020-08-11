function data = bb_MVPA_Independence_Test(data, subNum, activation)
% Independence can be proved at the IV(1) level.  Correlation between consecutive 
%  blocks require that the mean activation for each condition be subtracted.

% Determine when each condition was shown in the experiment sequence.
seq_unique = unique(data.sub(subNum).cond_sequence);
for cond = 1:length(seq_unique)
    condition(cond).vec = data.sub(subNum).cond_sequence == seq_unique(cond);
end            

% For each activation value, subtract the corresponding condition mean.  This is done within each voxel.
for cond = 1:length(condition)
    activation(condition(cond).vec,:) = activation(condition(cond).vec,:) - repmat(mean(mean(activation(condition(cond).vec,:))),[sum(condition(cond).vec),size(activation,2)]);
end

% Compute the correlation within each voxel
len = length(data.sub(subNum).cond_sequence);
corr = NaN(size(activation,2),1);
for jj=1:size(activation,2)
    R = corrcoef(activation(1:(len-1),jj),activation(2:len,jj));
    corr(jj) = R(1,2);
end
corr = corr(~isnan(corr));

data.sub(subNum).study.stats.iv(1).mvpa.searchlight.ind.median = median(corr);
data.sub(subNum).study.stats.iv(1).mvpa.searchlight.ind.corr_hist = hist(corr);
data.sub(subNum).study.stats.iv(1).mvpa.searchlight.ind.title  = 'Independence Test: Correlation Between Consecutive Blocks';
data.sub(subNum).study.stats.iv(1).mvpa.searchlight.ind.figfilename = sprintf('%s/%s_%s_Independence_Test.jpg', data.sub(subNum).qc.MVPAaaqcdir, data.sub(subNum).study.info.id, data.sub(subNum).info.name);
figure; data.sub(subNum).study.stats.iv(1).mvpa.searchlight.ind.corr_hist;
xlabel('Pairwise Block Correlation', 'FontSize', 24);
ylabel('Consecutive Pairs of Blocks', 'FontSize', 24);
saveas(figure(1), data.sub(subNum).study.stats.iv(1).mvpa.searchlight.ind.figfilename);
close all;
end