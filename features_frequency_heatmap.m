clear;
close;
%% load SVM results
m_dir = 'P:\Jose_Chonay\classification\SVM_ratios\';
task = 'mobi';
run = 'End';
n_ft = '15';
load([m_dir 'all-subs_' n_ft 'ft_' task '_' run '.mat']);
%% create matrix with labels
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
%% count frequency
% Initialize frequency matrix
freq_mat = zeros(size(matrix));

% Loop over each list in list_cell
for k = 1:size(SVM_results.features_lb,1)
    current_list = SVM_results.features_lb{k};
    
    % Loop over each string in the current list
    for i = 1:length(current_list)
        str = current_list{i};
        
        % Create a logical matrix where M equals str
        match_mat = strcmp(matrix, str);
        
        % Increment frequency matrix where match occurs
        freq_mat = freq_mat + match_mat;
    end
end
%% metrics
accuracy_means = mean(SVM_results.accuracy);
accuracy_sd    = std(SVM_results.accuracy);
AUC_means      = mean(SVM_results.AUC);
AUC_sd        = std(SVM_results.AUC);
pre_means     = mean(SVM_results.precission);
pre_sd        = std(SVM_results.precission);
re_means     = mean(SVM_results.recall);
re_sd          = std(SVM_results.recall);
%% plot
c_rgb = hex2rgb('#ea43b1');
ncolors = 256; % Number of colors in the map
cmap = [linspace(1, c_rgb(1), ncolors)', ...
        linspace(1, c_rgb(2), ncolors)', ...
        linspace(1, c_rgb(3), ncolors)'];

% Plot heatmap
figure;
imagesc(freq_mat);
colormap(cmap);
tick_pos = [3 8 13 18];

% Set X ticks
xticks(tick_pos);
xticklabels(ROI);
set(gca, 'FontSize', 20);  % Set font size for ticks
ax = gca;   % current axes
ax.XRuler.TickLabelGapOffset = 18;   % negative pushes labels further down
ax.YRuler.TickLabelGapOffset = 18;   % negative pushes labels further left

% Set Y ticks
yticks(tick_pos);
yticklabels(ROI);
set(gca, 'FontSize', 20);  % Set font size for ticks
ax = gca;   % current axes
ax.XRuler.TickLabelGapOffset = 20;   % negative pushes labels further down
ax.YRuler.TickLabelGapOffset = 20;   % negative pushes labels further left

c = colorbar;
clim([0 ceil(max(freq_mat(:)))]);
max_freq = max(freq_mat(:));
c.Ticks = 0:1:max_freq;

% add greek labels
greek_labels = {'\theta', '\alpha', '\beta', '\gamma', '\gamma'''};

% X-axis Greek labels BELOW the FM/PM/LT/RT
for i = 1:20
    idx = mod(i-1, length(greek_labels)) + 1;
    text(i,20.4,greek_labels{idx}, ...  % adjust "4.4" for spacing below RT row
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',16);
end

% Y-axis Greek labels LEFT of FM/PM/LT/RT
for i = 1:20
    idx = mod(i-1, length(greek_labels)) + 1;
    text(0.3,i,greek_labels{idx}, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontSize',16);
end
title(sprintf(['Frequency of features kept during RFE ( ' n_ft ' features set as max) for estimating an SVM classifier \n ' ...
    'for MTLR and control groups in task ' task ' during run ' run]));
%%
print(gcf,['P:\Jose_Chonay\classification\SVM_ratios\feature_frequency_' n_ft '_' task '_' run '.png'],'-dpng','-r1500');