%%
clear
close all;
%% load data
run = 'End';
task = 'mobi';
ROI = {'FM'; 'PM'; 'LT'; 'RT'};

%%
tables_ptnts = {};
tables_ctrls = {};
for i_r = 1:size(ROI,1)
    load(['P:\Sein_Jeung\Project_Watermaze\WM_EEG_Results\BOSC\BOSC_probe_' task '_' run '_' ROI{i_r} '.mat']);
    trials_patients = {};
    for i_p = 1:10
        trials_p = squeeze(mean(boscOutputs{1,i_p}.pepisode)); % average over electrode groups
        trials_patients{i_p} = trials_p; % trials X features
    end
    
    fBandsIdx = boscOutputs{1,1}.config.eBOSC.F;
    patients =   vertcat(trials_patients{:});

    bands_patients = util_WM_bandratios(patients, ROI{i_r}, fBandsIdx);
    tables_ptnts{i_r} = bands_patients;

    trials_controls = {};
    for i_c = 11:30
        trials_c = squeeze(mean(boscOutputs{1,i_c}.pepisode));
        trials_controls{i_c-10} = trials_c;
    end
    
    controls = vertcat(trials_controls{:});

    bands_controls = util_WM_bandratios(controls, ROI{i_r}, fBandsIdx);
    tables_ctrls{i_r} = bands_controls;
end
%% ratio matrix
patients_concat = horzcat(tables_ptnts{:});
controls_concat = horzcat(tables_ctrls{:});

patients_mat = meanLogRatio(patients_concat);
controls_mat = meanLogRatio(controls_concat);

ratio_mat = struct();
ratio_mat.patients = patients_mat;
ratio_mat.controls = controls_mat;

%% save
save_dir = 'P:\Jose_Chonay\classification\logit_ratios\';
save([save_dir 'all-subs_' task '_' run '.mat'], "ratio_mat");
%%
figure; 
clims = [-2, 2];
n = 64;
% Color 1: #ff314c → [1.0000, 0.1922, 0.2980]
c1 = [1.0000, 0.1922, 0.2980];  % red-pink
% White
white = [1, 1, 1];
% Color 2: #e68b00 → [0.9020, 0.5451, 0.0000]
c2 = [0.9020, 0.5451, 0.0000];  % orange
% Gradient from red-pink to white
half1 = [linspace(c1(1), white(1), n)', ...
         linspace(c1(2), white(2), n)', ...
         linspace(c1(3), white(3), n)'];
% Gradient from white to orange
half2 = [linspace(white(1), c2(1), n)', ...
         linspace(white(2), c2(2), n)', ...
         linspace(white(3), c2(3), n)'];

% Combine the two halves
cmap = [half1; half2];

% Apply the colormap
colormap(cmap);
subplot(1,2,1); imagesc(squeeze(mean(ratio_mat.patients,3)), clims); title('Mobile MTLR'); 
xticks([3:5:18]);  yticks([3:5:18]);  xticklabels(ROI); yticklabels(ROI); 
set(gca, 'FontSize', 20);  % Set font size for ticks
ax = gca;   % current axes
ax.XRuler.TickLabelGapOffset = 15;   % negative pushes labels further down
ax.YRuler.TickLabelGapOffset = 15;   % negative pushes labels further left


subplot(1,2,2); imagesc(squeeze(mean(ratio_mat.controls,3)), clims); title('Mobile CTRL'); 
xticks([3:5:18]);  yticks([3:5:18]);  xticklabels(ROI); yticklabels(ROI); 
set(gca, 'FontSize', 20);  % Set font size for ticks
ax = gca;   % current axes
ax.XRuler.TickLabelGapOffset = 15;   % negative pushes labels further down
ax.YRuler.TickLabelGapOffset = 15;   % negative pushes labels further left
colorbar;
% add greek labels
greek_labels = {'\theta', '\alpha', '\beta', '\gamma', '\gamma'''};

% X-axis Greek labels BELOW the FM/PM/LT/RT
for i = 1:20
    idx = mod(i-1, length(greek_labels)) + 1;
    subplot(1,2,1)
    text(i,20.4,greek_labels{idx}, ...  % adjust "4.4" for spacing below RT row
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',16);
    subplot(1,2,2)
        text(i,20.4,greek_labels{idx}, ...  % adjust "4.4" for spacing below RT row
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',16);
end

% Y-axis Greek labels LEFT of FM/PM/LT/RT
for i = 1:20
    idx = mod(i-1, length(greek_labels)) + 1;
    subplot(1,2,1)
    text(0.3,i,greek_labels{idx}, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontSize',16);
    subplot(1,2,2)
        text(0.3,i,greek_labels{idx}, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontSize',16);
end
%%
print(gcf,[save_dir 'all-subs_' task '_' run '.png'],'-dpng','-r1500');