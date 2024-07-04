function data = S_Searchlight_Analysis(data, subNum, activation, debug)
for i=1:2
slsize = i; %Searchlight's surrounding layer thickness
threshp = 0.0001;

% Determine the seachlight filename.  If it already exists, exit this function.
if ~debug
    data.sub(subNum).mvpadir = sprintf('%s/%s/mvpa_results', data.sub(subNum).study.info.studydir, data.sub(subNum).info.name);
    if (exist(data.sub(subNum).mvpadir, 'dir') ~= 7)
        mkdir(data.sub(subNum).mvpadir);
    end
    location = sprintf('%s/%s/mvpa_results/%s', data.sub(subNum).study.info.studydir, data.sub(subNum).info.name,data.sub(subNum).study.stats.iv(2).tag);
    if (exist(location, 'dir') ~= 7)
        mkdir(location);
    end
    temp = strrep(sprintf('%0.4f', threshp), '.', 'p');
    data.sub(subNum).study.stats.iv(2).mvpa.searchlight_file = sprintf('%s/%s_Searchlight_Task_Effect_anal%s_SLsize_%d.nii', location, data.sub(subNum).info.name, temp, slsize);
    if exist(data.sub(subNum).study.stats.iv(2).mvpa.searchlight_file, 'file') == 2
        return;
    end
end

%Use "len" to ensure that each class is equally represented.
len = min(size(activation.stim,4),size(activation.null,4));
len2 = len*2;
data.sub(subNum).study.stats.iv(2).mvpa.dataframe_len = len2;
% dataframe = [activation.stim;activation.null];
%dataframe = dataframe - repmat(mean(dataframe,1),len*2,1);
% dataframe = zscore(dataframe);

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
data.sub(subNum).study.stats.iv(2).mvpa.anal.threshv = binoinv((1-threshp),tst_len*2,0.5)/(tst_len*2);

% Run the searchlight
fprintf('%s  Starting Searchlight Calcs\n', datestr(now)); if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
for k=1:k_max
    for j=1:j_max
        for i=1:i_max
            % Check that the voxel activation values are non-zero
            if all(activation.stim(i,j,k,:)) && all(activation.null(i,j,k,:))
                dvoxel=[];
                for k_ind = max((k-slsize),1):min((k+slsize),k_max)
                    for j_ind = max((j-slsize),1):min((j+slsize),j_max)
                        for i_ind = max((i-slsize),1):min((i+slsize),i_max)
                            % Check that the voxel activation values are non-zero
                            if all(activation.stim(i_ind,j_ind,k_ind,:)) && all(activation.null(i_ind,j_ind,k_ind,:))
                                temp = [squeeze(activation.stim(i_ind,j_ind,k_ind,:)); squeeze(activation.null(i_ind,j_ind,k_ind,:))];
                                dvoxel = [dvoxel, temp];
                            end
                        end
                    end
                end
                %Perform the classification analysis using cross validation.
                correct = zeros(tst_len,1);
                for r = 1:(tst_len)
                    rr = tst(r);
                    trainvect = logical(ones(len2,1)); trainvect(rr) = 0; trainvect(len+rr) = 0;
                    tstvect = logical(zeros(len2,1)); tstvect(rr) = 1; tstvect(len+rr) = 1;
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
fprintf('%s  Finished Searchlight Calcs\n', datestr(now));if ~debug;diary('/autofs/cluster/vincent/programs/matlab/suraj/diary_Test_Part2');end;
% Threshhold the accuracy values
voxelacc(voxelacc < data.sub(subNum).study.stats.iv(2).mvpa.anal.threshv) = 0;
% Check whether there was any significant result
if ~all(voxelacc == 0)
    data.sub(subNum).study.stats.iv(2).mvpa.anal.sig = 1;
end
% Assign the 98% accuracy to threshmaxv
nz = sort(nonzeros(voxelacc));
data.sub(subNum).study.stats.iv(2).mvpa.anal.threshmaxv = nz(ceil(length(nz)*0.98));
% Read in template file and modify it to create the searchlight file
if debug
    data.sub(subNum).study.stats.iv(2).mvpa.voxelacc = voxelacc;
else
    template = sprintf('%s/%s/atlas/%s_fmr_fmrspace_ave.nii', ...
        data.sub(subNum).study.info.studydir, data.sub(subNum).info.name, data.sub(subNum).info.name);
    mri = MRIread(template);
    mri.vol = voxelacc;
    mri.fspec = data.sub(subNum).study.stats.iv(2).mvpa.searchlight_file;
    MRIwrite(mri, mri.fspec);
    clear mri;
end
clear voxelacc;
end
end
