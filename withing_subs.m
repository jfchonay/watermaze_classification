clear;
clc;
%%
% config structure used for selection ROI
config_param.chanGroups(1).key           = 'FM';
config_param.chanGroups(1).full_name     = 'Frontal-midline';
config_param.chanGroups(1).chan_names    = {'y1','y2','y3','y25','y32'}; 
config_param.chanGroups(2).key           = 'PM';
config_param.chanGroups(2).full_name     = 'Parietal-midline';
config_param.chanGroups(2).chan_names    = {'r9', 'r10', 'r11', 'r27', 'r32'}; 
config_param.chanGroups(3).key           = 'LT';
config_param.chanGroups(3).full_name     = 'Left-temporal';
config_param.chanGroups(3).chan_names    = {'g1', 'y16', 'r15', 'r13'}; 
config_param.chanGroups(4).key           = 'RT';
config_param.chanGroups(4).full_name     = 'Right-temporal';
config_param.chanGroups(4).chan_names    = {'g24','y20', 'r18', 'r20'}; 

%% Load results structs
band = 'theta';
run = 'end';
main_dir = ['P:\Jose_Chonay\frequency_sliding\' band];
% load populations and task
pt_s = load([main_dir '\ptnts_probe_stat_fsliding_' run '.mat']);
ct_s = load([main_dir '\ctrls_probe_stat_fsliding_' run '.mat']);

pt_m = load([main_dir '\ptnts_probe_mobi_fsliding_' run '.mat']);
ct_m = load([main_dir '\ctrls_probe_mobi_fsliding_' run '.mat']);

%% cut 200 ms from the start and end
to_cut = 0.2*ct_s.ctrl_struct.srate; % in samples
time = ct_s.ctrl_struct.ori_time{1,1}(1+to_cut:end-to_cut);

patients_s = pt_s.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_s = ct_s.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

patients_m = pt_m.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_m = ct_m.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

% cells of the data nested in condition then by population
control = {controls_s, controls_m};
patient = {patients_s, patients_m};
trials = {control, patient};

permr = struct();
permr.params.dependent_samples = 1;
permr.params.p_threshold = 0.05;
permr.params.two_sided = 1;

%%
for pop = 1:2 % control or patients
    for Fi = 1:4
        chan_indices = find(cellfun(@(x) any(strcmp(x, config_param.chanGroups(Fi).chan_names)), ct_m.ctrl_struct.chan_labels));
        switch pop
            case 1 % controls
                % run a permutation within subjects
                x = squeeze(mean(trials{1,pop}{1,1}(chan_indices, :, :),1)); % stationary
                y = squeeze(mean(trials{1,pop}{1,2}(chan_indices, :, :),1)); % mobile
                % paramaters set to within subject analysis
                [clusters, p_values, t_sums] = permutest(x,y,permr.params.dependent_samples,permr.params.p_threshold,[],permr.params.two_sided);
                %save results to structure
                permr.controls(Fi).ROI_key = config_param.chanGroups(Fi).key;
                permr.controls(Fi).ROI_full_name = config_param.chanGroups(Fi).full_name;
                permr.controls(Fi).clusters = clusters;
                permr.controls(Fi).p_values = p_values;
                permr.controls(Fi).t_sums = t_sums;
                permr.controls(Fi).clusters_astime = cellfun(@(idx) time(idx), clusters, 'UniformOutput', false);
            case 2 % patients
                % run a permutation test between subjects
                x = squeeze(mean(trials{1,pop}{1,1}(chan_indices, :, :),1)); % patients
                y = squeeze(mean(trials{1,pop}{1,2}(chan_indices, :, :),1)); % control
                % paramaters set to between subject analysis and p
                % threshold of 0.05 and two sided test
                [clusters, p_values, t_sums] = permutest(x,y,permr.params.dependent_samples,permr.params.p_threshold,[],permr.params.two_sided);
                % save results to structure
                permr.patients(Fi).ROI_key = config_param.chanGroups(Fi).key;
                permr.patients(Fi).ROI_full_name = config_param.chanGroups(Fi).full_name;
                permr.patients(Fi).clusters = clusters;
                permr.patients(Fi).p_values = p_values;
                permr.patients(Fi).t_sums = t_sums;
                permr.patients(Fi).clusters_astime = cellfun(@(idx) time(idx), clusters, 'UniformOutput', false);

        end
    end
end
%% save the structure
save(fullfile(main_dir, ['within_subs_perm_' run '.mat']), "permr");