%%
clear
close all;
%% load data
mDir = 'P:\Jose_Chonay\classification\logit_ratios\';
task = 'mobi';
run = 'End';
%%
load([mDir 'all-subs_' task '_' run '.mat']);
%%
ROI = {'FM'; 'PM'; 'LT'; 'RT'};
bands = {'theta','alpha', 'beta', 'gamma', 'highgamma'};
% Create all 20 combinations
combinations = cell(20, 1);
index = 1;
for letter = 1:4
    for number = 1:5
        combinations{index} = [ROI{letter} bands{number}];
        index = index + 1;
    end
end

% Create the 20×20 combination matrix
matrix = cell(20, 20);
for row = 1:20
    row_base = combinations{row}; 
    for col = 1:20
        col_base = combinations{col}; 
        % Combine row base with column base
        matrix{row, col} = [row_base '_' col_base]; 
    end
end

%%
patients = permute(ratio_mat.patients, [3 2 1]);
controls = permute(ratio_mat.controls, [3 2 1]);
nFeaturesToKeep = 15;
k = 5; % Number of folds for cross-validation
SVM_results = struct();
%%
%for all 190 features
%% for all triangle features
% Logical index for lower triangle without diagonal
mask = tril(true(size(patients,2)), -1);  % lower triangle, excluding diagonal
num_elements = nnz(mask);  % number of elements in the lower triangle without diagonal

% Preallocate result matrix: each row is one frame, each column one lower-tri element
X_patients = zeros(size(patients,1), num_elements);
X_controls = zeros(size(controls,1), num_elements);

% Loop over each slice
for i_p = 1:size(patients,1)
    mat_p = squeeze(patients(i_p,:,:));     % Get the i-th 20x20 matrix
    X_patients(i_p,:) = mat_p(mask);  % Extract lower triangle (no diagonal)
end
for i_c = 1:size(controls,1)
    mat_c = squeeze(controls(i_c,:,:));
    X_controls(i_c,:) = mat_c(mask);
end

Y_patients = repmat((1), 240, 1); % patients = 1
Y_controls = repmat((-1), 476, 1); % controls = -1
ft_labels = matrix(mask)';
%%
trials = [X_controls; X_patients];
labels = [Y_controls; Y_patients];

rng(1); % For reproducibility

% Preallocate performance metrics
models_all = cell(k,1);
features_all = cell(k,1);
accuracy_all = zeros(k, 1);
F1_all = zeros(k, 1);
precision_all = zeros(k, 1);
recall_all = zeros(k, 1);
AUC_all = zeros(k, 1);
ft_labels_all = cell(k,1);

cv = cvpartition(labels, 'KFold', k);

%% 
if isempty(gcp('nocreate'))
    parpool;  % Uses default pool
end

%%
parfor fold = 1:k
    fprintf('\nFold %d/%d\n', fold, k);

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

    %% Feature selection using backward sequential selection

    % Loss function for SVM
    fun = @(XTrain, yTrain, XTest, yTest) ...
        loss(fitcsvm(XTrain, yTrain, ...
        'KernelFunction', 'linear', ...
        'Weights', get_weights(yTrain)), ...
        XTest, yTest);

   opts = statset('Display','off','UseParallel',false); % Turn off display for parallel speed

    % Inner cross-validation partition (use holdout to avoid nested KFold parallel)
    cv_inner = cvpartition(Y_train, 'KFold', 10);

    % Disable parallel within sequentialfs
    [inmodel, ~] = sequentialfs(fun, X_train, Y_train, ...
        'cv', cv_inner, ...
        'direction', 'backward', ...
        'nfeature', nFeaturesToKeep,...
        'options', opts); % Important for parfor compatibility

    selectedFeatures = find(inmodel);

    %% Train final model on selected features
    finalModel = fitcsvm(X_train(:, selectedFeatures), Y_train, ...
        'KernelFunction', 'linear', ...
        'Weights', get_weights(Y_train), ...
        'Standardize', true, ...
        'OptimizeHyperparameters', 'auto', ...
        'HyperparameterOptimizationOptions', ...
        struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
               'UseParallel', false)); % Disable nested parallel

    %% Predict and evaluate
    [Y_pred, scores] = predict(finalModel, X_test(:, selectedFeatures));

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
    accuracy_all(fold) = accuracy;
    F1_all(fold) = F1;
    precision_all(fold) = precision;
    recall_all(fold) = recall;
    AUC_all(fold) = AUC;
    features_all{fold,1} = selectedFeatures;
    ft_labels_all{fold,1} = ft_labels(selectedFeatures);
    models_all{fold} = finalModel;
end

%% repeated features
% Concatenate all arrays (works with different sizes)
allNumbers = horzcat(features_all{:});

uniqueNumbers = unique(allNumbers);
counts = zeros(size(uniqueNumbers));

% Count how many arrays each number appears in
for i_u = 1:length(uniqueNumbers)
    num = uniqueNumbers(i_u);
    % Count occurrences in each cell using cellfun
    counts(i_u) = sum(cellfun(@(x) ismember(num, x), features_all));
end

% Find numbers appearing in ≥threshold arrays
threshold = 3;  % You can change this as needed
repeatingNumbers = uniqueNumbers(counts >= threshold);
repeatingFeatures = ft_labels(repeatingNumbers);
%% save relevant data
SVM_results.modality    = task;
SVM_results.run         = run;
SVM_results.accuracy    = accuracy_all;
SVM_results.precission  = precision_all;
SVM_results.recall      = recall_all;
SVM_results.F1          = F1_all;
SVM_results.AUC         = AUC_all;
SVM_results.ft_used     = nFeaturesToKeep;
SVM_results.features_lb = ft_labels_all;
SVM_results.repeat_ft   = repeatingFeatures;
SVM_results.models      = models_all;
%%
save_dir = ['P:\Jose_Chonay\classification\SVM_ratios\' 'all-subs_' num2str(nFeaturesToKeep) 'ft_' task '_' run '.mat']; 
save(save_dir, "SVM_results");