%%
clear
close all;
%% load data
mDir = 'P:\Jose_Chonay\classification\logit_ratios\';
task = 'mobi';
run = 'Start';
%%
load([mDir 'all-subs_' task '_' run '.mat']);
patients = permute(ratio_mat.patients, [3 2 1]);
controls = permute(ratio_mat.controls, [3 2 1]);

%% for 19 features

% ROI = {'FM'; 'PM'; 'LT'; 'RT'};
% band = {'theta','alpha', 'beta', 'gamma', 'highgamma'};
% 
% idx = (find(strcmp('RT', ROI))*5) - (5 - find(strcmp('theta', band)));
% 
% for i_p = 1:size(patients,1)
%     f_p = squeeze(patients(i_p,idx,:))';
%     f_p(idx) = []; 
%     p_ft(i_p,:) = f_p;
% end
% 
% for i_c = 1:size(controls,1)
%     f_c = squeeze(controls(i_c, idx,:))';
%     f_c(idx) = [];
%     c_ft(i_c,:) = f_c;
% end
% 
% X_patients = p_ft;
% X_controls = c_ft;
% 
% Y_patients = repmat((1), 240, 1); % patients = 1
% Y_controls = repmat((-1), 476, 1); % controls = -1
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
%% Prepare data

trials = [X_controls; X_patients];
labels = [Y_controls; Y_patients];

rng(1); % For reproducibility

nFeaturesToKeep = 26;
k = 5; % Number of folds for cross-validation
cv = cvpartition(labels, 'KFold', k);
% Initialize parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool; % Uses default number of workers
end

% Preallocate results structure
ft_cv = struct('accuracy', cell(1, length(nFeaturesToKeep:-2:1)), ...
               'AUC', [], 'pre', [], 're', []);
%%
%% 
parfor nf_idx = 1:length(nFeaturesToKeep:-2:1)
    nf = nFeaturesToKeep - (nf_idx-1)*2;  % Calculate actual nFeatures value
    
    % Initialize temporary variables for this iteration
    temp_accuracy = zeros(1, k);
    temp_F1 = zeros(1, k);
    temp_precision = zeros(1, k);
    temp_recall = zeros(1, k);
    temp_AUC = zeros(1, k);
    temp_features = cell(k, 1);
    
    for fold = 1:k
        fprintf('\nProcessing nFeatures=%d, Fold %d/%d\n', nf, fold, k);
        
        trainIdx = training(cv, fold);
        testIdx = test(cv, fold);
        
        X_train = trials(trainIdx, :);
        Y_train = labels(trainIdx);
        X_test = trials(testIdx, :);
        Y_test = labels(testIdx);
        
        %% Compute class weights
        class1_count = sum(Y_train == 1);
        class2_count = sum(Y_train == -1);
        weight_for_1 = 1 / class1_count;
        weight_for_minus1 = 1 / class2_count;
        get_weights = @(y) (y == 1)*weight_for_1 + (y == -1)*weight_for_minus1;
        
        %% Feature selection
        fun = @(XTrain, yTrain, XTest, yTest) ...
            loss(fitcsvm(XTrain, yTrain, ...
            'KernelFunction', 'linear', ...
            'Weights', get_weights(yTrain)), ...
            XTest, yTest);
        
        opts = statset('Display', 'final');
        cv_inner = cvpartition(Y_train, 'KFold', 5);
        
        [inmodel, ~] = sequentialfs(fun, X_train, Y_train, ...
            'cv', cv_inner, ...
            'direction', 'backward', ...
            'nfeatures', nf, ...
            'options', opts);
        
        selectedFeatures = find(inmodel);
        
        %% Train final model
        finalModel = fitcsvm(X_train(:, selectedFeatures), Y_train, ...
            'KernelFunction', 'linear', ...
            'Weights', get_weights(Y_train));
        
        %% Predict and evaluate
        [Y_pred, scores] = predict(finalModel, X_test(:, selectedFeatures));
        
        [~, ~, ~, AUC] = perfcurve(Y_test, scores(:, 2), 1);
        
        TP = sum((Y_test == 1) & (Y_pred == 1));
        FP = sum((Y_test == -1) & (Y_pred == 1));
        FN = sum((Y_test == 1) & (Y_pred == -1));
        
        precision = TP / (TP + FP + eps);
        recall = TP / (TP + FN + eps);
        F1 = 2 * (precision * recall) / (precision + recall + eps);
        accuracy = sum(Y_pred == Y_test) / numel(Y_test);
        
        % Store metrics in temporary variables
        temp_accuracy(fold) = accuracy;
        temp_F1(fold) = F1;
        temp_precision(fold) = precision;
        temp_recall(fold) = recall;
        temp_AUC(fold) = AUC;
        temp_features{fold} = selectedFeatures;
    end
    
    % Assign results to output structure
    ft_cv(nf_idx).accuracy = temp_accuracy;
    ft_cv(nf_idx).AUC = temp_AUC;
    ft_cv(nf_idx).pre = temp_precision;
    ft_cv(nf_idx).re = temp_recall;
    ft_cv(nf_idx).features = temp_features;
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
threshold = 5;  % You can change this as needed
repeatingNumbers = uniqueNumbers(counts >= threshold);
%% get means and std
accuracy_means = zeros(size(ft_cv,2),1);
accuracy_sd    = zeros(size(ft_cv,2),1);
AUC_means      = zeros(size(ft_cv,2),1);
AUC_sd         = zeros(size(ft_cv,2),1);
pre_means      = zeros(size(ft_cv,2),1);
pre_sd         = zeros(size(ft_cv,2),1);
re_means       = zeros(size(ft_cv,2),1);
re_sd          = zeros(size(ft_cv,2),1);
n_features     = nFeaturesToKeep:-2:1;
for n_s = 1:size(ft_cv,2)
    accuracy_means(n_s) = mean(ft_cv(n_s).accuracy);
    accuracy_sd(n_s)    = std(ft_cv(n_s).accuracy);
    AUC_means(n_s)      = mean(ft_cv(n_s).AUC);
    AUC_sd(n_s)         = std(ft_cv(n_s).AUC);
    pre_means(n_s)      = mean(ft_cv(n_s).pre);
    pre_sd(n_s)         = std(ft_cv(n_s).pre);
    re_means(n_s)       = mean(ft_cv(n_s).re);
    re_sd(n_s)          = std(ft_cv(n_s).re);
end
%% plot
% Create figure
figure;
hold on;

% Plot Accuracy with error bars
errorbar(n_features, accuracy_means, accuracy_sd, ...
    'o--', 'LineWidth', 1.5, 'MarkerSize', 8, 'CapSize', 10, ...
    'DisplayName', 'Accuracy');

% Plot AUC with error bars
errorbar(n_features, AUC_means, AUC_sd, ...
    's--', 'LineWidth', 1.5, 'MarkerSize', 8, 'CapSize', 10, ...
    'DisplayName', 'AUC');

errorbar(n_features, pre_means, pre_sd, ...
    'o--', 'LineWidth', 1.5, 'MarkerSize', 8, 'CapSize', 10, ...
    'DisplayName', 'Precission');

errorbar(n_features, re_means, re_sd, ...
    's--', 'LineWidth', 1.5, 'MarkerSize', 8, 'CapSize', 10, ...
    'DisplayName', 'Recall');

% Customize the plot
set(gca, 'XDir', 'reverse');  % Reverse x-axis (since we go from 19 →2 features)
xlabel('Number of Features');
ylabel('Metric Value');
title(sprintf(['Evaluation metrics by number of features used for classifying controls and MTLR \n' ...
    ' for the start of mobile task']));
legend('show');
grid on;
%%
print(gcf,'P:\Jose_Chonay\classification\feature_selection\.png','-dpng','-r1500');