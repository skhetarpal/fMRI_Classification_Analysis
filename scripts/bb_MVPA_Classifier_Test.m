function results = bb_MVPA_Classifier_Test(group1, group2, classifier)

%Use "len" to ensure that each class is equally represented.
results.len = min(size(group1,1),size(group2,1));
len = results.len;
if size(group1,1) ~= size(group2,1)
    perm1 = randperm(size(group1,1));
    perm2 = randperm(size(group2,1));
    sorted1 = sort(perm1(1:len));
    sorted2 = sort(perm2(1:len));
    group1 = group1(sorted1,:);
    group2 = group2(sorted2,:);
end
d=[group1; group2];
for jj=1:size(d,2)
    d(:,jj) = zscore(d(:,jj));
end

%Build the classifier. Use Cross-Validation to determine accuracy.
runvect = 1:len*2;
condvect = [ones(len, 1); ones(len, 1)*2]; %These are the machine learning labels
results.cond1_acc=zeros(1,len); results.cond2_acc=zeros(1,len);
fprintf('%s  Determining ROI accuracy\n', datestr(now));
for rr=1:len
    trainvect = runvect~=rr&runvect~=(len+rr);
    tstvect = runvect==rr|runvect==(len+rr);
    trainmat = d(trainvect,:);
    tstmat = d(tstvect,:);
    traincond = condvect(trainvect);
    tstcond = condvect(tstvect);
    if strcmp(classifier, 'svm')
        SVM = svmtrain(trainmat,traincond);
        predictions = svmclassify(SVM, tstmat);
    end
    if strcmp(classifier, 'classify')
        predictions = classify(tstmat, trainmat, traincond, 'diaglinear');
    end
    results.cond1_acc(rr) = predictions(1)==tstcond(1);
    results.cond2_acc(rr) = predictions(2)==tstcond(2);
end
results.cond1_mean_acc = mean(results.cond1_acc);
results.cond2_mean_acc = mean(results.cond2_acc);
results.acc = mean([results.cond1_mean_acc, results.cond2_mean_acc]);
results.acc_stderr = std([results.cond1_acc results.cond2_acc]) / sqrt(2*len);

% Permutation Test to Determine P-Value.
perms = 100;
permutation_acc1 = zeros(1,perms);
permutation_acc2 = zeros(1,perms);
fprintf('%s  Determining ROI permutation accuracies\n', datestr(now));
for pp=1:perms
    correct1=zeros(len,1);
    correct2=zeros(len,1);
    condvectprime = condvect(randperm(numel(condvect))); %Scramble the labels
    for rr=1:len
        trainvect = runvect~=rr&runvect~=(len+rr);
        tstvect = runvect==rr|runvect==(len+rr);
        trainmat = d(trainvect,:);
        tstmat = d(tstvect,:);
        traincond = condvectprime(trainvect);
        tstcond = condvectprime(tstvect);
        if strcmp(classifier, 'svm')
            SVM = svmtrain(trainmat,traincond);
            predictions = svmclassify(SVM, tstmat);
        end
        if strcmp(classifier, 'classify')
            predictions = classify(tstmat, trainmat, traincond, 'diaglinear');
        end
        correct1(rr) = predictions(1)==tstcond(1);
        correct2(rr) = predictions(2)==tstcond(2);
    end
    permutation_acc1(pp) = mean(correct1);
    permutation_acc2(pp) = mean(correct2);
end
fprintf('%s  Finished ROI permutation accuracies\n', datestr(now));
permutation_acc = sort(mean([permutation_acc1; permutation_acc2]));
results.pval_cond1 = (sum(permutation_acc1 > results.cond1_mean_acc)+1)/(perms+1);
results.pval_cond2 = (sum(permutation_acc2 > results.cond2_mean_acc)+1)/(perms+1);
results.pval = (sum(permutation_acc > results.acc)+1)/(perms+1);
results.perm_mean = mean(permutation_acc);
results.perm_std = std(permutation_acc);
results.perm_stderr = results.perm_std / sqrt(perms);
results.perm_upper = permutation_acc(round(perms*.95));
results.perm_lower = permutation_acc(round(perms*.05));
end
