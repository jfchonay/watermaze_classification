%%
clear
close all;
%% load data
task = 'mobi';
run = 'Start';
ROI = 'RT';
load(['P:\Sein_Jeung\Project_Watermaze\WM_EEG_Results\BOSC\BOSC_probe_' task '_' run '_' ROI '.mat']);

%% separate into trials of patients and controls
% define frequency band
s_band = 4;
e_band = 8;
patients_cells = {};
for i_p = 1:10
    trials_p = squeeze(mean(boscOutputs{1,i_p}.pepisode,1)); % average over electrode groups
    % find indices of the frequency bands we want to extract
    idx_s = find(boscOutputs{1,i_p}.config.eBOSC.F >= s_band, 1);
    idx_e = find(boscOutputs{1,i_p}.config.eBOSC.F <= e_band, 1,'last');
    patients_cells{i_p}= trials_p(:,idx_s:idx_e); % trials X features
end
patients_trials = vertcat(patients_cells{:});
% repeat procedure for controls
controls_cells = {};
for i_c = 11:30
    trials_c = squeeze(mean(boscOutputs{1,i_c}.pepisode,1));
    % find indices of the frequency bands we want to extract
    idx_sc = find(boscOutputs{1,i_c}.config.eBOSC.F >= s_band, 1);
    idx_ec = find(boscOutputs{1,i_c}.config.eBOSC.F <= e_band, 1,'last');
    controls_cells{i_c-10} = trials_c(:,idx_sc:idx_ec);
end
controls_trials = vertcat(controls_cells{:});

Y_patients = repmat((1), 240, 1); % patients = 1
Y_controls = repmat((-1), 476, 1); % controls = -1

%% prepare data for SVM
trials = [mean(controls_trials,2); mean(patients_trials,2)];
labels = [Y_controls; Y_patients];

%% k fold partition
cv = cvpartition(size(trials,1), 'KFold', 5);

% metric storage
allAccuracy = zeros(cv.NumTestSets,1);
allPrecision = zeros(cv.NumTestSets,1);
allRecall = zeros(cv.NumTestSets,1);
allF1 = zeros(cv.NumTestSets,1);
allAUC = zeros(cv.NumTestSets,1);
%% 
if isempty(gcp('nocreate'))
    parpool;  % Uses default pool
end
%%
parfor fold = 1:cv.NumTestSets
    fprintf('\nFold %d/%d\n', fold, cv.NumTestSets);

    trainIdx = training(cv, fold);
    testIdx = test(cv, fold);

    X_train = trials(trainIdx, :);
    Y_train = labels(trainIdx);
    X_test = trials(testIdx, :);
    Y_test = labels(testIdx);

    % Compute class weights from training set
    class1_count = sum(Y_train == 1);
    class2_count = sum(Y_train == -1);
    weight_for_1 = 1 / class1_count;
    weight_for_minus1 = 1 / class2_count;
    get_weights = @(y) (y == 1)*weight_for_1 + (y == -1)*weight_for_minus1;

    %% Train final model on selected features
    finalModel = fitcsvm(X_train, Y_train, ...
        'KernelFunction', 'linear', ...
        'Weights', get_weights(Y_train), ...
        'Standardize', true, ...
        'OptimizeHyperparameters', 'auto', ...
        'HyperparameterOptimizationOptions', ...
        struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
               'UseParallel', false)); % Disable nested parallel

    %% Predict and evaluate
    [Y_pred, scores] = predict(finalModel, X_test);

    % Use score for class 1 for ROC
    [~, ~, ~, AUC] = perfcurve(Y_test, scores(:, 2), 1);

    % Metrics
    TP = sum((Y_test == 1) & (Y_pred == 1));
    FP = sum((Y_test == -1) & (Y_pred == 1));
    FN = sum((Y_test == 1) & (Y_pred == -1));

    precision = TP / (TP + FP + eps);
    recall = TP / (TP + FN + eps);
    F1 = 2 * (precision * recall) / (precision + recall + eps);
    accuracy = sum(Y_pred == Y_test) / numel(Y_test);

    % Store metrics
    allAccuracy(fold) = accuracy;
    allF1(fold) = F1;
    allPrecision(fold) = precision;
    allRecall(fold) = recall;
    allAUC(fold) = AUC;
end
%% visualize it
labelsm = string(cell2mat(labels(:,1)));
decisionBoundary = -SVMModel.Bias / SVMModel.Beta;


figure;
hold on;
scatter(trials(labelsm == 'controls'), zeros(sum(labelsm == 'controls'), 1), 'b', 'filled');
scatter(trials(labelsm == 'patients'), zeros(sum(labelsm == 'patients'), 1), 'r', 'filled');
line([decisionBoundary, decisionBoundary], [-0.1, 0.1], 'Color', 'g', 'LineWidth', 2);

% Highlight the support vectors
scatter(SVMModel.SupportVectors, zeros(size(SVMModel.SupportVectors, 1), 1), 100, 'k', 'x', 'LineWidth', 0.5);

% Add labels and legend
xlabel('Feature');
legend('controls', 'patients','Decision Boundary', 'Support Vectors');
title('SVM with 1 Feature: Support Vectors');
hold off;
%% save results
SVMresults = struct();
SVMresults.model = SVMModel;
SVMresults.data.trials = trials;
SVMresults.data.labels = labels;
SVMresults.data.low_band = s_band;
SVMresults.data.upper_band = e_band;
SVMresults.cv = cv;
SVMresults.metrics.precision = allPrecision;
SVMresults.metrics.recall = allRecall;
SVMresults.metrics.f1 = allF1;
SVMresults.metrics.auc = allAUC;

to_save = ['P:\Jose_Chonay\classification\SVM_eBOSC\SVMresults_theta_' task '_' run '_' ROI '.mat'];
save(to_save, "SVMresults");