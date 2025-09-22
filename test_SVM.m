clear;
clc;
%% Load results 
run = 'Start';
main_dir = 'P:\Jose_Chonay\classification\power_results';

% load populations and task
power_table = readtable([main_dir '\resultsTable_power_probe_'  run '.xlsx']);

% break table into groups and conditions
rf = rowfilter(power_table);
controls_stat = power_table(rf.Group == '2' & rf.Setup == '1',:);
patients_stat = power_table(rf.Group == '1' & rf.Setup == '1',:);

controls_mobi = power_table(rf.Group == '2' & rf.Setup == '2',:);
patients_mobi = power_table(rf.Group == '1' & rf.Setup == '2',:);

%% load trials
% band = 'alpha';
% 
% main_dir = ['P:\Jose_Chonay\frequency_sliding\' band];
% 
% session = 'probe';
% task = 'mobi';
% event = 'fsliding';
% ROI = 'RT';
% run = 'start';
% 
% controls_mobi = load([main_dir '\ctrls_' session '_' task '_' event '_' ROI '_' run '.mat']);
% patients_mobi = load([main_dir '\ptnts_' session '_' task '_' event '_' ROI '_' run '.mat']);
%% organizing labels and data
labels = [repmat({'controls'}, length(controls_mobi.Group), 1) ; repmat({'patients'}, length(patients_mobi.Group), 1)];
theta = [controls_mobi.RTtheta ; patients_mobi.RTtheta];
gamma = [controls_mobi.FMgamma ; patients_mobi.FMgamma];
trials = [theta, gamma];
%DTW = [controls_mobi.DTW ; patients_mobi.DTW];
%features = [power, DTW];

% %z_data = zscore(power, [], 1);

%% organize labels and data
% labels = [controls_mobi.trialstrct.labels(:,1);  patients_mobi.trialstrct.labels(:,1)];
% trials = [controls_mobi.trialstrct.trials; patients_mobi.trialstrct.trials];

%% cv hold out
cv = cvpartition(length(labels), 'HoldOut', 0.2);
Xtrain = trials(training(cv),:);
Ytrain = labels(training(cv));
Xtest = trials(test(cv),:);
Ytest = labels(test(cv));

%% set weights
% we want to balance the weights to give the minority class a more balanced
% approach
weights = calculateWeights(Ytrain);

%% cv kfold
% 
% cv = cvpartition(size(trials,1), 'KFold', 5);
% 
% % metric storage
% allPrecision = zeros(cv.NumTestSets,1);
% allRecall = zeros(cv.NumTestSets,1);
% allF1 = zeros(cv.NumTestSets,1);
% allLoss = zeros(cv.NumTestSets,1);


%% perform k fold
% for i = 1:cv.NumTestSets
%     trainIdx = training(cv, i); % 
%     testIdx = test(cv, i); %
% 
%     Xtrain = trials(trainIdx, :);
%     Ytrain = labels(trainIdx);
%     Xtest = trials(testIdx, :);
%     Ytest = labels(testIdx);
% 
%     % we want to balance the weights to give the minority class a more balanced
%     % approach
%     weights = calculateWeights(Ytrain);
% 
%     % Train SVM Model for This Fold
%     SVMModel = fitcsvm(Xtrain, Ytrain, ...
%                        'KernelFunction', 'linear', ...
%                        'Weights', weights,...
%                        'Standardize', true, ...
%                        'OptimizeHyperparameters', 'auto', ...
%                        'HyperparameterOptimizationOptions', ...
%                        struct('AcquisitionFunctionName', 'expected-improvement-plus',...
%                        'MaxObjectiveEvaluations', 50));
% 
%     % Make Predictions for This Fold
%     Ypred = predict(SVMModel, Xtest);
% 
%     % Compute Confusion Matrix Values
%     TP = sum((cell2mat(Ytest) == 'patients') & (cell2mat(Ypred) == 'patients'));
%     FP = sum((cell2mat(Ytest) == 'controls') & (cell2mat(Ypred) == 'patients'));
%     FN = sum((cell2mat(Ytest) == 'patients') & (cell2mat(Ypred) == 'controls'));
% 
%     % Compute Performance Metrics for This Fold
%     Precision = TP / (TP + FP); % Positive predictive value, when a patient is correctly classified
%     Recall = TP / (TP + FN); % True positive rate, few false negatives
%     F1 = 2 * (Precision * Recall) / (Precision + Recall); % balance of precission and recall
%     Loss = loss(SVMModel, Xtest, Ytest); 
% 
%     % Store Results
%     allPrecision(i) = Precision;
%     allRecall(i) = Recall;
%     allF1(i) = F1;
%     allLoss(i) = Loss;
% end

%% model fitting hold out
SVMModel = fitcsvm(Xtrain, Ytrain, ...
                   'KernelFunction', 'linear', ...
                   'Weights', weights,...
                   'Standardize', true, ...
                   'OptimizeHyperparameters', 'auto', ...
                   'HyperparameterOptimizationOptions', ...
                   struct('AcquisitionFunctionName', 'expected-improvement-plus',...
                   'MaxObjectiveEvaluations', 50));

%% evaluation on cv
CVModel = crossval(SVMModel, 'Leaveout', 'on');
cvLoss = kfoldLoss(CVModel);
%% predicted Y and scores
[Ypred, scores] = predict(SVMModel, Xtest);

%% 
% Confusion Matrix
figure;
confusionchart(Ytest, Ypred);
title('Confusion Matrix');

% Compute AUC-ROC
[Xroc, Yroc, T, AUC] = perfcurve(Ytest, scores(:,2), 'controls');
fprintf('AUC-ROC: %.3f\n', AUC);

% Plot ROC Curve
figure;
plot(Xroc, Yroc, 'b', 'LineWidth', 2);
xlabel('False Positive Rate'); ylabel('True Positive Rate');
title('ROC Curve');
grid on;

%% metrics
% Precision, Recall, F1-score
TP = sum((cell2mat(Ytest) == 'patients') & (cell2mat(Ypred) == 'patients'));
FP = sum((cell2mat(Ytest) == 'controls') & (cell2mat(Ypred) == 'patients'));
FN = sum((cell2mat(Ytest) == 'patients') & (cell2mat(Ypred) == 'controls'));
Precision = TP / (TP + FP);
Recall = TP / (TP + FN);
F1 = 2 * (Precision * Recall) / (Precision + Recall);

fprintf('Precision: %.3f\n', Precision);
fprintf('Recall: %.3f\n', Recall);
fprintf('F1-score: %.3f\n', F1);

