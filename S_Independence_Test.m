function data = S_Independence_Test(data, subNum, activation)
% Independence can be proved at the IV(2) level using correlation between consecutive task effect values.

image_size = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims;
% Compute the correlation within each voxel
len = size(activation.stim,4);
corr = NaN(data.sub(subNum).study.mristudy.fmri.mvpa.image_dims);
for k=1:image_size(3)
    for j=1:image_size(2)
        for i=1:image_size(1)
            R = corrcoef(activation.stim(i,j,k,1:(len-1)),activation.stim(i,j,k,2:len));
            corr(i,j,k) = R(1,2);
        end
    end
end
corr = reshape(corr, image_size(1)*image_size(2)*image_size(3), 1);
corr = corr(~isnan(corr));

data.sub(subNum).study.stats.iv(2).mvpa.searchlight.ind.median = median(corr);
figure; hist(corr);
data.sub(subNum).study.stats.iv(2).mvpa.searchlight.ind.title  = 'Independence Test: Correlation Between Consecutive Blocks';
data.sub(subNum).study.stats.iv(2).mvpa.searchlight.ind.figfilename = ...
    sprintf('%s/%s_%s_Independence_Test.jpg', data.sub(subNum).qc.MVPAaaqcdir, ...
    data.sub(subNum).study.info.id, data.sub(subNum).info.name);
xlabel('Pairwise Block Correlation', 'FontSize', 24);
ylabel('Consecutive Pairs of Blocks', 'FontSize', 24);
saveas(gcf, data.sub(subNum).study.stats.iv(2).mvpa.searchlight.ind.figfilename);
close all;
end