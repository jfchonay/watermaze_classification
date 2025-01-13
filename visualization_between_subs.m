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

pt_s = load([main_dir '\ptnts_probe_stat_fsliding_' run '.mat']);
ct_s = load([main_dir '\ctrls_probe_stat_fsliding_' run '.mat']);

pt_m = load([main_dir '\ptnts_probe_mobi_fsliding_' run '.mat']);
ct_m = load([main_dir '\ctrls_probe_mobi_fsliding_' run '.mat']);

%% cut data

to_cut = 0.2*ct_s.ctrl_struct.srate;
time = ct_s.ctrl_struct.ori_time{1,1}(1+to_cut:end-to_cut);

patients_s = pt_s.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_s = ct_s.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

patients_m = pt_m.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_m = ct_m.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

stat = {patients_s, controls_s};
mobi = {patients_m, controls_m};
trials = {stat, mobi};
colors = {'#C75C22', '#D74674', '#A762BB', '#1781D7'};

for c = 1:size(colors,2)
    colors_rgb{c} = hex2rgb(colors{c});
end

%% create the figure

figure;
for cond = 1:2 % stat and mobi
    for Fi = 1:4 % ROI
        switch cond
            case 1
                subplot(4, 2, Fi*2 - 1);
                
                % find channel indices
                chan_indices = find(cellfun(@(x) any(strcmp(x, config_param.chanGroups(Fi).chan_names)), ct_m.ctrl_struct.chan_labels));
                
                for Gi = 1:2 % group ind 1 = patients 2 = controls
                    % Extract data for the channel of interest (e.g., channel Fi)
                    channel_data = squeeze(mean(trials{1,cond}{1,Gi}(chan_indices, :, :),1)); % [1151×10 or X20 double]
                    
                    % Compute mean and confidence interval (95%) across
                    % participants for one band
                    mean_data = mean(channel_data, 2); % [1151×1]
                    std_data = std(channel_data, 0, 2); % [1151×1]
                    n_participants = size(channel_data, 2);
                    ci = 1.96 * (std_data / sqrt(n_participants)); % Confidence interval
                    
                    % Plot the mean with confidence interval as a ribbon
                    hold on;
                    fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{Gi}, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence interval ribbon
                    plot(time, mean_data, 'Color', colors_rgb{Gi}, 'LineWidth', 1.5); % Mean
                    
                    % Labels and formatting
                    title(['Channel ' config_param.chanGroups(Fi).key ]);
                    %legend('','patient','','control');
                    ylim([5 7]);
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                end
                % run a permutation test to get some error bars
                x = squeeze(mean(trials{1,cond}{1,1}(chan_indices, :, :),1));
                y = squeeze(mean(trials{1,cond}{1,2}(chan_indices, :, :),1));
                [clusters, p_values, t_sums, permutation_distribution ] = permutest(x,y,0,0.05);
                if ~isempty(t_sums)
                    for i_c = 1:length(clusters)
                        hold on
                        y_value = 5.2; % Set the y-coordinate
                        plot(time(clusters{1,i_c}), y_value * ones(size(time(clusters{1,i_c}))), 'k.', 'LineWidth', 0.2); % 
                    end
                end
            case 2
                 subplot(4, 2, Fi*2);
                
                % find channel indices
                chan_indices = find(cellfun(@(x) any(strcmp(x, config_param.chanGroups(Fi).chan_names)), ct_m.ctrl_struct.chan_labels));
                
                for Gi = 1:2 % group ind 1 = patients 2 = controls
                    % Extract data for the channel of interest (e.g., channel Fi)
                    channel_data = squeeze(mean(trials{1,cond}{1,Gi}(chan_indices, :, :),1)); % [1251×10 double]
                    
                    % Compute mean and confidence interval (95%) across
                    % participants for one band
                    mean_data = mean(channel_data, 2); % [1251×1]
                    std_data = std(channel_data, 0, 2); % [1251×1]
                    n_participants = size(channel_data, 2);
                    ci = 1.96 * (std_data / sqrt(n_participants)); % Confidence interval
                    
                    % Plot the mean with confidence interval as a ribbon
                    hold on;
                    fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{Gi+2}, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence interval ribbon
                    plot(time, mean_data, 'Color', colors_rgb{Gi+2}, 'LineWidth', 1.5); % Mean
                    
                    % Labels and formatting
                    title(['Channel ' config_param.chanGroups(Fi).key ]);
                    %legend('','patient','','control');
                    ylim([5 7]);
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                end
                 % run a permutation test to get some error bars
                x = squeeze(mean(trials{1,cond}{1,1}(chan_indices, :, :),1));
                y = squeeze(mean(trials{1,cond}{1,2}(chan_indices, :, :),1));
                [clusters, p_values, t_sums, permutation_distribution ] = permutest(x,y,0,0.05);
                if ~isempty(t_sums)
                    for i_c = 1:length(clusters)
                        hold on
                        y_value = 5.2; % Set the y-coordinate
                        plot(time(clusters{1,i_c}), y_value * ones(size(time(clusters{1,i_c}))), 'k.', 'LineWidth', 0.2); % 
                    end
                end

            sgtitle(sprintf(['Instantenous frequency for all ' run ' trials filtered for ' ...
                band ' band \n static vs mobile']));
        end
    end
end
%print(gcf,[main_dir '\between_subs_' run '.png'],'-dpng','-r1500');