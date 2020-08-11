function datavar = bb_MVPA_Searchlight_Analysis(data, subNum, tag, datavar, threshp, name, matrix1, matrix2)

slsize = 2; %Searchlight's surrounding layer thickness

% Determine the seachlight filename.  If it already exists, exit this function.
temp = sprintf('%0.4f', threshp);
temp = strrep(temp, '.', 'p');
sl_name = sprintf('%s_Searchlight_%s_anal%s_SLsize_%d', data.sub(subNum).info.name, name, temp, slsize);
datavar.searchlight_file = sprintf('%s/%s/mvpa_results/%s/%s', data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, tag, sl_name);
fspec = sprintf('%s.nii', datavar.searchlight_file);
if exist(fspec, 'file') == 2
    return;
end

%% Compute Searchlight Classification Accuracy

%Use "len" to ensure that each class is equally represented.
len = min(size(matrix1,1),size(matrix2,1));
len2 = len*2;
datavar.dataframe_len = len2;
if size(matrix1,1) ~= size(matrix2,1)
%     matrix1 = matrix1(randperm(len),:);
%     matrix2 = matrix2(randperm(len),:);
    matrix1 = matrix1(1:len,:);
    matrix2 = matrix2(1:len,:);
end
dataframe = [matrix1;matrix2];
%dataframe = dataframe - repmat(mean(dataframe,1),len*2,1);
dataframe = zscore(dataframe);

condvect = [ones(len, 1); ones(len, 1)*2]; %These are the machine learning labels
i_max = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(1);
j_max = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(2);
k_max = data.sub(subNum).study.mristudy.fmri.mvpa.image_dims(3);
voxelacc = zeros(i_max, j_max, k_max);

%For the sake of speed, limit the number of blocks that will be tested for accuracy.
fl = floor(len/10);
if fl~=0
    tst_len = 10;
    tst = fl:fl:(fl*10);
elseif fl == 0
    tst_len = len;
    tst = 1:len;
end

%Determine the minimum accuracy theshhold
datavar.anal.threshv = binoinv((1-threshp),tst_len*2,0.5)/(tst_len*2);

% Run the searchlight
fprintf('%s  Starting Searchlight Calcs\n', datestr(now)); diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
for k=1:1:k_max
    for j=1:1:j_max
        for i=1:1:i_max
            voxel = i + (j-1)*i_max + (k-1)*i_max*j_max;
            if all(dataframe(:,voxel))
                % If the voxel index contains data, it is in the brain, so create an matrix containing voxels within the searchlight.
                dvoxel=[];
                for k_ind = max((k-slsize),1):min((k+slsize),k_max)
                    for j_ind = max((j-slsize),1):min((j+slsize),j_max)
                        for i_ind = max((i-slsize),1):min((i+slsize),i_max)
                            i_voxel = i_ind + (j_ind-1)*i_max + (k_ind-1)*i_max*j_max;
                            if all(dataframe(:,i_voxel))
                                dvoxel = [dvoxel, [squeeze(dataframe(:,i_voxel))]];
                            end
                        end
                    end
                end
                %Perform the classification analysis using cross validation.
                correct = zeros(tst_len,1);
                for r = 1:(tst_len)
                    rr = tst(r);
                    trainvect = ones(len2); trainvect(rr) = 0; trainvect(len+rr) = 0;
                    tstvect = zeros(len2); tstvect(rr) = 1; tstvect(len+rr) = 1;
                    trainmat = dvoxel(trainvect,:);
                    tstmat = dvoxel(tstvect,:);
                    traincond = condvect(trainvect);
                    tstcond = condvect(tstvect);
                    class = classify(tstmat, trainmat, traincond, 'diaglinear');
                    correct(r) = mean(class==tstcond);
%                     SVM = svmtrain(trainmat,traincond);
%                     predictions = svmclassify(SVM, tstmat);
%                     correct(rr) = mean(predictions==tstcond);
                end
                voxelacc(i,j,k) = mean(correct); %Store the mean accuracy
            end
        end
    end
end
fprintf('%s  Finishing Searchlight Calcs\n', datestr(now));diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');
% Threshhold the accuracy values
voxelacc(voxelacc < datavar.anal.threshv) = 0;
% Check whether there was any significant result
if ~all(voxelacc == 0)
    datavar.anal.sig = 1;
end
% Assign the 98% accuracy to threshmaxv
nz = sort(nonzeros(voxelacc));
datavar.anal.threshmaxv = nz(ceil(length(nz)*0.98));
% Read in template file and modify it to create the searchlight file
template = sprintf('%s/%s/spm_results/IV2_Task/%s_spmT_Task_UNC0p001_neg.nii', ...
    data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
mri = MRIread(template);
mri.vol = voxelacc;
clear voxelacc;
mri.fspec = fspec;
MRIwrite(mri, mri.fspec);
clear mri;
datavar = bb_MVPA_Searchlight_Atlas_Transformation(data, subNum, datavar);
datavar = bb_MVPA_Create_Searchlight_Overlays(data, subNum, datavar, threshp);
end